%% Analyze_all_responses.m
% Aggregates response data across experiments and performs population analysis.
%
% Prerequisites:
%   - Run Find_EPSPs.m for each experiment first
%   - amp_zscores.mat and fake_amp_zscores.mat in intermediate/ (auto-generated
%     there on first run if not already present)

% --- Configuration ---
paths = jsondecode(fileread('config/paths.json'));
data_dir = paths.mosaix_data_dir;
output_dir = fullfile(pwd, paths.output_dir);

% Experiment groupings
thalamus = {'0730', '0802'};
contra = {'0811', '0814'};

% Stimulus timing windows
PreStim = 1:32;
Stim = 34:39;
PostStim = 43:63;
PostEPSP = 175:200;
trials = 1:100;

%% Load and aggregate all response files from experiment output folders

all_experiments = [thalamus, contra];
response_table = table();
for i = 1:numel(all_experiments)
    exp_id = all_experiments{i};
    fn = fullfile(output_dir, exp_id, sprintf('Responses_by_type_%s.mat', exp_id));
    if ~exist(fn, 'file')
        warning('Response file not found for experiment %s, skipping.', exp_id);
        continue;
    end
    response_data = load(fn);
    response_table = vertcat(response_table, response_data.to_write);
end
fprintf('Loaded %d experiments, %d total cells\n', numel(all_experiments), size(response_table, 1));

%% Compute response amplitudes

response_amps = cellfun(@(x) mean(x(1, PostStim, trials), [2 3]), response_table.snip_data, 'UniformOutput', false);
response_amps = cell2mat(response_amps);
norm_amps = response_amps ./ response_table.cheriff_data;

%% Z-score based responder detection

% Load or compute bootstrap z-scores for successive trial averaging
intermediate_dir = fullfile(pwd, paths.intermediate_dir);
if ~exist(intermediate_dir, 'dir')
    mkdir(intermediate_dir);
end

amp_zscores_path = fullfile(intermediate_dir, 'amp_zscores.mat');
fake_amp_zscores_path = fullfile(intermediate_dir, 'fake_amp_zscores.mat');

if exist(amp_zscores_path, 'file') && exist(fake_amp_zscores_path, 'file')
    fprintf('Loading pre-computed z-scores from intermediate/\n');
    amp_zscores = load(amp_zscores_path);
    amp_zscores = amp_zscores.amp_zscores;
    fake_amp_zscores = load(fake_amp_zscores_path);
    fake_amp_zscores = fake_amp_zscores.fake_amp_zscores;
else
    fprintf('Computing z-scores of successively averaged trials (this may take a while)...\n');
    sampleSpace = PreStim;
    n_trials_range = 20:190;

    amp_zscores = zeros(numel(n_trials_range), size(response_table, 1));
    boot_std = zeros(numel(n_trials_range), size(response_table, 1));
    boot_means = zeros(numel(n_trials_range), size(response_table, 1));
    resp_amps = zeros(numel(n_trials_range), size(response_table, 1));
    fake_amp_zscores = zeros(numel(n_trials_range), 500, size(response_table, 1));

    for i = 1:size(response_table, 1)
        traces = response_table.snip_data{i};
        for j = 1:numel(n_trials_range)

            mean_trace = mean(traces(:,:,1:19+j), 3);
            mean_baseline = mean(mean_trace(sampleSpace));
            boot_means(j,i) = mean_baseline;
            std_baseline = std(mean_trace(sampleSpace));
            boot_std(j,i) = std_baseline/sqrt(numel(PostStim));
            resp_amp = mean(mean_trace(PostStim));
            resp_amps(j,i) = resp_amp;
            amp_z = (resp_amp - boot_means(j,i))/boot_std(j,i);
            amp_zscores(j,i) = amp_z;

            for k = 1:500
                rand_fake_idx = randsample(sampleSpace, numel(PostStim), true);
                mean_fake_resp = mean(mean_trace(rand_fake_idx));
                fake_z = (mean_fake_resp - boot_means(j,i))/boot_std(j,i);
                fake_amp_zscores(j,k,i) = fake_z;
            end
        end
    end

    save(amp_zscores_path, 'amp_zscores');
    save(fake_amp_zscores_path, 'fake_amp_zscores');
    fprintf('Saved z-scores to intermediate/\n');
end

is_responding = any((amp_zscores <= -3.75), 1);
fake_is_responding = any((fake_amp_zscores <= -3.75), 1);
fake_pos_rate = sum(fake_is_responding, 2) / 500;
mean_fpr = mean(fake_pos_rate);

responders_table = response_table(logical(is_responding), :);
non_responders_table = response_table(~logical(is_responding), :);

%% Determine optimal trial count per cell (z-score convergence)

first_maxes = zeros(sum(is_responding), 1);
norm_amp_zs = amp_zscores(:, logical(is_responding));
responder_idx = find(is_responding);

for i = 1:numel(responder_idx)
    z_scores_to_use = medfilt1(-amp_zscores(:, responder_idx(i)));
    max_min = (max(z_scores_to_use) - min(z_scores_to_use));
    norm_amp_zs(:, i) = (z_scores_to_use - min(z_scores_to_use)) ./ max_min;
    indices = find(norm_amp_zs(:, i) > 0.85);
    first_maxes(i) = quantile(indices, 0.9);
end

cell_max_zscores = zeros(numel(first_maxes), 1);
for i = 1:numel(first_maxes)
    cell_max_zscores(i) = amp_zscores(floor(first_maxes(i)), responder_idx(i));
end

[sorted, idx] = sort(first_maxes, 'ascend');
cm = [0 0 0; 0.8 0.5 0.5; 1 1 1];
cmi = interp1([0; 0.8; 1], cm, (0:0.0001:1));
figure; imagesc(medfilt1(norm_amp_zs(:, idx)')); colormap(cmi); colorbar; hold on;
plot(first_maxes(idx), 1:numel(first_maxes), 'b.');
title('Z-score convergence by cell');

%% Trim trials to optimal count per cell

zscore_index = first_maxes;
for i = 1:size(responders_table, 1)
    trials_data = responders_table.snip_data{i};
    end_index = floor(zscore_index(i));
    select_trials = trials_data(:, :, 1:end_index);
    responders_table.snip_data(i) = {select_trials};
end

% Filter out direct CheRiff-positive cells
responders_table = responders_table(responders_table.cheriff_pos == 0, :);

%% Recompute amplitudes for filtered responders

response_amps = cellfun(@(x) mean(x(1, PostStim, :), [2 3]), responders_table.snip_data, 'UniformOutput', false);
response_amps = cell2mat(response_amps);
norm_amps = response_amps ./ responders_table.cheriff_data;

% Median normalize responders (by experiment and slice)
exps = unique(responders_table.experiment);
med_norm_amps = zeros(size(response_amps));
med_amps = zeros(numel(exps), 1);

for i = 1:numel(exps)
    in_exp = strcmp(responders_table.experiment, exps{i});
    for j = unique(responders_table.slice_data(logical(in_exp)))'
        sub_amps = response_amps(logical(in_exp & responders_table.slice_data == j));
        med_amp = median(sub_amps);
        med_amps(i) = med_amp;
        mm_amp = sub_amps ./ med_amp;
        med_norm_amps(logical(in_exp & responders_table.slice_data == j)) = mm_amp;
    end
end

%% Build combined table (responders + non-responders)

responders_table.("is_responder") = ones(size(responders_table, 1), 1);
non_responders_table.("is_responder") = zeros(size(non_responders_table, 1), 1);
all_table = vertcat(responders_table, non_responders_table);

all_response_amps = cellfun(@(x) mean(x(1, PostStim, :), [2 3]), all_table.snip_data, 'UniformOutput', false);
all_response_amps = cell2mat(all_response_amps);
all_norm_amps = all_response_amps ./ all_table.cheriff_data;

% Median normalize all cells to the responder median
all_exps = unique(all_table.experiment);
all_med_norm_amps = zeros(size(all_response_amps));
for i = 1:numel(all_exps)
    in_exp = strcmp(all_table.experiment, all_exps{i});
    for j = unique(all_table.slice_data(logical(in_exp)))'
        sub_amps = all_response_amps(logical(in_exp & all_table.slice_data == j));
        mm_amp = sub_amps ./ med_amps(i);
        all_med_norm_amps(logical(in_exp & all_table.slice_data == j)) = mm_amp;
    end
end

%% Save output tables

save(fullfile(output_dir, 'all_table.mat'), 'all_table', 'all_med_norm_amps');
save(fullfile(output_dir, 'responders_table.mat'), 'responders_table');
save(fullfile(output_dir, 'non_responders_table.mat'), 'non_responders_table');

fprintf('Analysis complete.\n');
fprintf('  Total cells: %d\n', size(response_table, 1));
fprintf('  Responders: %d\n', size(responders_table, 1));
fprintf('  Non-responders: %d\n', size(non_responders_table, 1));
fprintf('  Mean false positive rate: %.4f\n', mean_fpr);
fprintf('  Saved tables to %s\n', output_dir);
