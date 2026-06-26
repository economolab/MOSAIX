cd '/Volumes/SynMap Data';

file_list = dir(pwd);
file_names = string({file_list.name}');

response_re = regexp(file_names, '^Response');
response_idx = find(~cellfun(@isempty,response_re));

fn = fullfile(pwd, file_names(response_idx(1),:));
response_data = load(fn);
response_table = response_data.to_write;
for i = 2:length(response_idx)
    fn = fullfile(pwd, file_names(response_idx(i),:));
    response_data = load(fn);
    response_table = vertcat(response_table, response_data.to_write);
end

%response_table = response_table(response_table.cheriff_pos == 0, :);
%L1_filter = response_table.L1_cheriff > 0.7;
%response_table = response_table(logical(L1_filter), :);

thalamus = {'0730', '0802'};
contra = {'0811', '0814'}; %'1008'

PreStim = 1:32;
Stim = 34:39;
PostStim = 43:63;
PostEPSP = 175:200;
%LongPostEPSP = 130:180;
trials = 1:100; %Full series should be something like 180...

response_amps = cellfun(@(x) mean(x(1,PostStim,trials), [2 3]), response_table.snip_data, 'UniformOutput', false);
response_amps = cell2mat(response_amps);
norm_amps = response_amps./response_table.cheriff_data; % Normalized to cheriff columnar value

%% SNR requirement determination if you use a set number of trials %%

% mean_trace = cellfun(@(x) mean(x(:,:,trials), 3), response_table.snip_data, 'UniformOutput', false);
% mean_trace = cell2mat(mean_trace);
% 
% sampleSpace = PreStim;
% 
% amp_zscores = zeros(size(response_table, 1), 1);
% boot_std = zeros(size(response_table, 1), 1);
% boot_means = zeros(500, size(response_table, 1));
% 
% for i = 1:size(response_table, 1)
%     for j = 1:500
%         rand_idx = randsample(sampleSpace, numel(PostStim));
%         mean_baseline = mean(mean_trace(i,rand_idx));
%         boot_means(j,i) = mean_baseline;
%     end
%     boot_std(i) = std(boot_means(:,i));
%     amp_z = (norm_amps(i) - mean(boot_means(:,i)))/boot_std(i);
%     amp_zscores(i) = amp_z;
% end
% 
% %orig_std = std(mean_trace(:,PreStim), 0, 2)/sqrt(21);
% %figure; scatter(orig_std, boot_std);
% %figure; histogram(boot_means(:,300), 100); xline(detection_lims(300));
% 
% detection_lims = mean(boot_means,1)' - (2.576 * boot_std); %2.576 = 0.99%, 3.719 = 0.9999%
% figure; histogram(detection_lims, 50);
% 
% figure;
% scatter(detection_lims, norm_amps, 65, 'filled'); hold on;
% plot([-10 0], [-10 0]);
% ylim([-0.5 0.2]);
% xlim([-0.5 0]);
% ylabel("'Response' amplitude");
% xlabel('Detection limit');
% xline(-0.2, '--', {'response = ~0.75mV'});
% %yline(-0.1, '--', {'response = ~0.4mV'});
% %yline(-2.576, '--', {'Z = -2.576'});
% %yline(-1.96, '--', {'Z = -1.0'});
% 
% sub_zscores = amp_zscores(find(amp_zscores >= 2.576));
% sub_std = boot_std(find(amp_zscores >= 2.576));
% % figure;
% % scatter(sub_std, sub_zscores, 65, 'filled');
% % %xline(-0.25, '--', {'response = ~1mV'});
% % yline(2.576, '--', {'Z = 2.576'});

%to_keep = detection_lims >= -0.2; % -0.13 is ~0.5 mV, -0.25 is about 1mV
%response_table = response_table(to_keep,:);

%% Select responders based on detection limit cutoff

% to_keep = detection_lims >= -0.2; % -0.13 is ~0.5 mV, -0.25 is about 1mV
% cells_responding = logical((norm_amps <= detection_lims) & to_keep);
% cells_not_responding = logical((norm_amps > detection_lims) & to_keep);
% 
% responders_table = response_table(logical(cells_responding & ~response_table.cheriff_pos),:);
% non_responders_table = response_table(logical(cells_not_responding & ~response_table.cheriff_pos),:);
% response_table = response_table(to_keep,:);
% 
% response_amps = cellfun(@(x) mean(x(1,PostStim,trials), [2 3]), responders_table.snip_data, 'UniformOutput', false);
% response_amps = cell2mat(response_amps);
% norm_amps = response_amps./responders_table.cheriff_data; % Normalized to cheriff columnar value
% 
% % Min-max normalize (by slice)
% exps = unique(responders_table.experiment);
% mm_norm_amps = zeros(size(norm_amps));
% for i = 1:numel(exps)
%     in_exp = strcmp(responders_table.experiment, exps{i});
%     for j = unique(responders_table.slice_data(logical(in_exp)))'
%         sub_amps = norm_amps(logical(in_exp & responders_table.slice_data == j));
%         max_amp = -max(-sub_amps);
%         min_amp = -min(-sub_amps);
%         med_amp = median(sub_amps);
%         mm_amp = sub_amps/med_amp;
%         %mm_amp = (sub_amps-min_amp)/(max_amp-min_amp);
%         mm_norm_amps(logical(in_exp & responders_table.slice_data == j)) = mm_amp;
%     end
% end

%% SNR requirement determination if you don't know how many trials to use %%

% This block is commented because z_score data are saved. Uncomment to
% re-run
% sampleSpace = PreStim;
% 
% amp_zscores = zeros(numel(20:190), size(response_table, 1), 1);
% boot_std = zeros(numel(20:190), size(response_table, 1), 1);
% boot_means = zeros(numel(20:190), size(response_table, 1));
% resp_amps = zeros(numel(20:190), size(response_table, 1));
% fake_amp_zscores = zeros(numel(20:190), 500, size(response_table, 1), 1);
% 
% for i = 1:size(response_table, 1)
%     traces = response_table.snip_data{i};
%     for j = 1:numel(20:190)
% 
%         mean_trace = mean(traces(:,:,1:19+j), 3);
%         mean_baseline = mean(mean_trace(sampleSpace));
%         boot_means(j,i) = mean_baseline;
%         std_baseline = std(mean_trace(sampleSpace));
%         boot_std(j,i) = std_baseline/sqrt(numel(PostStim));
%         resp_amp = mean(mean_trace(PostStim));
%         resp_amps(j,i) = resp_amp;
%         amp_z = (resp_amp - boot_means(j,i))/boot_std(j,i);
%         amp_zscores(j,i) = amp_z;
% 
%         for k = 1:500
%             rand_fake_idx = randsample(sampleSpace, numel(PostStim), true);
%             mean_fake_resp = mean(mean_trace(rand_fake_idx));
%             fake_z = (mean_fake_resp - boot_means(j,i))/boot_std(j,i);
%             fake_amp_zscores(j,k,i) = fake_z;
%         end
%     end
% end

amp_zscores = load('amp_zscores.mat');
fake_amp_zscores = load('fake_amp_zscores.mat');
amp_zscores = amp_zscores.amp_zscores;
fake_amp_zscores = fake_amp_zscores.fake_amp_zscores;

is_responding = any((amp_zscores <= -3.75),1);
fake_is_responding = any((fake_amp_zscores <= -3.75),1);
fake_pos_rate = sum(fake_is_responding,2)/500;

responders_table = response_table(logical(is_responding),:);
non_responders_table = response_table(~logical(is_responding),:);
%responders_std = boot_std(:,logical(is_responding));
%non_responders_std = boot_std(:,~logical(is_responding));

first_maxes = zeros(sum(is_responding),1);
norm_amp_zs = amp_zscores(:,logical(is_responding));
responder_idx = find(is_responding);
for i = 1:numel(responder_idx)
    % if amp_zscores(1,i) < 0
    %     z_scores_to_use = -amp_zscores(:,i);
    % else
    %     z_scores_to_use = amp_zscores(:,i);
    % end
    z_scores_to_use = medfilt1(-amp_zscores(:,responder_idx(i)));
    max_min = (max(z_scores_to_use) - min(z_scores_to_use));
    norm_amp_zs(:,i) = (z_scores_to_use - min(z_scores_to_use))./max_min;
    indices = find(norm_amp_zs(:,i) > 0.85);
    %first_maxes(i) = median(indices);
    first_maxes(i) = quantile(indices, 0.9);
end

cell_max_zscores = zeros(numel(first_maxes), 1);
for i = 1:numel(first_maxes)
    cell_max_zscores(i) = amp_zscores(floor(first_maxes(i)),responder_idx(i));
end
%figure; histogram(cell_max_zscores,100); xlim([-15 0.5]);
[sorted, idx] = sort(first_maxes, 'ascend');
cm = [0 0 0; 0.8 0.5 0.5; 1 1 1];
cmi = interp1([0; 0.8; 1], cm, (0:0.0001:1));
figure; imagesc(medfilt1(norm_amp_zs(:,idx)')); colormap(cmi); colorbar; hold on;
plot(first_maxes(idx), 1:numel(first_maxes), 'b.');

zscore_index = first_maxes;

for i = 1:size(responders_table,1)
    trials = responders_table.snip_data{i};
    end_index = zscore_index(i);
    select_trials = trials(:,:,1:end_index);
    responders_table.snip_data(i) = {select_trials};
end

%responder_std_coords = horzcat(floor(first_maxes), (1:numel(first_maxes))');
%responder_std_coords = sub2ind(size(responders_std), responder_std_coords(:,1), responder_std_coords(:,2));
%best_std = vertcat(responders_std(responder_std_coords), non_responders_std(171,:)');
%figure; histogram(best_std*-3.75, 40);

responders_table = responders_table(responders_table.cheriff_pos == 0,:);

response_amps = cellfun(@(x) mean(x(1,PostStim,:), [2 3]), responders_table.snip_data, 'UniformOutput', false);
response_amps = cell2mat(response_amps);
norm_amps = response_amps./responders_table.cheriff_data; % Normalized to cheriff columnar value

% Median normalize responders (by slice)
exps = unique(responders_table.experiment);
mm_norm_amps = zeros(size(response_amps));
med_amps = size(numel(numel(exps)),1);
for i = 1:numel(exps)
    in_exp = strcmp(responders_table.experiment, exps{i});
    for j = unique(responders_table.slice_data(logical(in_exp)))'
        sub_amps = response_amps(logical(in_exp & responders_table.slice_data == j));
        max_amp = -max(-sub_amps);
        min_amp = -min(-sub_amps);
        med_amp = median(sub_amps);
        med_amps(i) = med_amp;
        mm_amp = sub_amps./med_amp;
        %mm_amp = (sub_amps-min_amp)/(max_amp-min_amp);
        mm_norm_amps(logical(in_exp & responders_table.slice_data == j)) = mm_amp;
    end
end

responders_table.("is_responder") = ones(size(responders_table,1),1);
non_responders_table.("is_responder") = zeros(size(non_responders_table,1),1);
all_table = vertcat(responders_table, non_responders_table);
all_response_amps = cellfun(@(x) mean(x(1,PostStim,:), [2 3]), all_table.snip_data, 'UniformOutput', false);
all_response_amps = cell2mat(all_response_amps);
all_norm_amps = all_response_amps./all_table.cheriff_data; % Normalized to cheriff columnar value

% Median normalize to the responder median (by slice)
all_exps = unique(all_table.experiment);
all_mm_norm_amps = zeros(size(all_response_amps));
for i = 1:numel(all_exps)
    in_exp = strcmp(all_table.experiment, all_exps{i});
    for j = unique(all_table.slice_data(logical(in_exp)))'
        sub_amps = all_response_amps(logical(in_exp & all_table.slice_data == j));
        %med_amp = median(sub_amps);
        mm_amp = sub_amps./med_amps(i);
        %mm_amp = (sub_amps-min_amp)/(max_amp-min_amp);
        all_mm_norm_amps(logical(in_exp & all_table.slice_data == j)) = mm_amp;
    end
end

%% Plot the z-score of successively added stim to see attenuation: %%

uniqueTypes = unique(responders_table.classif_data, 'stable');
[uniqueTypes, idx] = sort(uniqueTypes);
uniqueCols = unique(responders_table.classif_cols, 'stable');
uniqueCols = uniqueCols(idx);
uniqueSlices = unique(responders_table.slice_data, 'sorted');

sampleSpace = [PreStim PostEPSP];
amp_zscores = zeros(size(uniqueTypes,1), 180);
matchExp = ismember(responders_table.experiment, contra);
all_first_maxes = zeros(sum(matchExp), 1);
position = 1;

figure;

for i = 1:size(uniqueSlices,1)
    %ax = subplot(ceil(length(uniqueTypes)/5), 5, i);
    ax = subplot(1, length(uniqueSlices), i);

    matchSlice = ismember(responders_table.slice_data, uniqueSlices(i));
    %matchType = strcmpi(responders_table.classif_data, uniqueTypes{i});
    pos_cells = find(matchExp & matchSlice);
    z_table = responders_table(pos_cells,:);
    
    amp_zs = zeros(size(z_table,1),numel(1:180));

    for j = 1:size(z_table, 1)
        cell_resp = squeeze(z_table.snip_data{j});
        for k = 1:180
            resp = mean(cell_resp(:,1:k), 2);
            resp_amp = mean(resp(PostStim));
            base_mean = mean(resp(sampleSpace));
            base_std = std(resp(sampleSpace), 0);
            amp_z = (resp_amp - base_mean);%/base_std;
            amp_zs(j,k) = amp_z;
        end
    end
    amp_zscores(i,:) = median(amp_zs, 1);
    first_maxes = zeros(size(amp_zs, 1), 1);
    for j = 1:size(amp_zs, 1)
        amp_zs(j,:) = (amp_zs(j,:) - min(amp_zs(j,:)))./(max(amp_zs(j,:))-min(amp_zs(j,:)));
        indices = find(amp_zs(j,:) < 0.15);
        first_maxes(j) = median(indices);
        all_first_maxes(pos_cells(j)) = median(indices);
        position = position + 1;
    end
    [sorted, idx] = sort(first_maxes, 'ascend');
    imagesc(-amp_zs(idx, :));
    xline([100])
    %plot(amp_zscores(i,:), 'color', hex2rgb(uniqueCols{i}), 'LineWidth', 2);
    %title(strcat(uniqueTypes{i}, ', n=', num2str(size(z_table,1))));
    title(strcat('Slice ', num2str(uniqueSlices(i)), ', n=', num2str(size(z_table,1))));
end

figure;
sampleSpace = [PreStim PostEPSP];
exp_array = {contra thalamus};
amp_zscores = zeros(2, 180);
position = 1;
cm = [0 0 0; 0.8 0.5 0.5; 1 1 1];
cmi = interp1([0; 0.8; 1], cm, (0:0.0001:1));

for i = 1:2
    ax = subplot(1, 2, i);
    
    matchExp = ismember(responders_table.experiment, exp_array{i});
    all_first_maxes = zeros(sum(matchExp), 1);
    pos_cells = find(matchExp);
    z_table = responders_table(pos_cells,:);
    
    amp_zs = zeros(size(z_table,1),numel(1:180));

    for j = 1:size(z_table, 1)
        cell_resp = squeeze(z_table.snip_data{j});
        for k = 1:180
            resp = -mean(cell_resp(:,1:k), 2);
            resp_amp = mean(resp(PostStim));
            base_mean = mean(resp(sampleSpace));
            base_std = std(resp(sampleSpace), 0);
            amp_z = (resp_amp - base_mean)/base_std;%
            amp_zs(j,k) = amp_z;
        end
    end
    amp_zscores(i,:) = median(amp_zs, 1);
    first_maxes = zeros(size(amp_zs, 1), 1);
    norm_amp_zs = amp_zs;
    for j = 1:size(amp_zs, 1)
        norm_amp_zs(j,:) = (amp_zs(j,:) - min(amp_zs(j,:)))./(max(amp_zs(j,:))-min(amp_zs(j,:)));
        indices = find(norm_amp_zs(j,:) > 0.85);
        first_maxes(j) = median(indices);
        all_first_maxes(pos_cells(j)) = median(indices);
        position = position + 1;
    end
    [sorted, idx] = sort(first_maxes, 'ascend');
    imagesc(norm_amp_zs(idx, :));
    %colormap(cmi);
    %colorbar;
    xline([100])
end

responder_std = boot_std(cells_responding & ~response_table.cheriff_pos);
figure; scatter(responder_std(matchExp), all_first_maxes, 'filled');

%% Response amplitudes plotted by laminar depth %%

plotting_table = all_table(logical(ismember(all_table.experiment, contra)'),:);
figure; % This will show response amps by depth, normalized for mean column cheriff
res_plot = scatter(log2(all_mm_norm_amps(logical(ismember(all_table.experiment, contra)'))), plotting_table.depth_data, 65, hex2rgb(plotting_table.classif_cols), 'filled'); 
hold on;
ylim([-1.0 0]); 
xlim([-0.5 8]);
xline([1.0]);
title(strcat('n =  ', num2str(size(plotting_table,1)), ' cells'))
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Exp', plotting_table.experiment);
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Slice', plotting_table.slice_data);
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Roi Num', plotting_table.roi_data);
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('CellType', plotting_table.classif_data);


%% Response amps plotted by cell type or other metadata %%

figure;
boxchart(categorical(responders_table.classif_data), mm_norm_amps)%, ...
    %'GroupByColor', responders_table.experiment);

% By slice number:
figure;
boxchart(categorical(round(responders_table.slice_data,2)), mm_norm_amps, ...
    'GroupByColor', responders_table.experiment);

% By cheriff
type_oi_pos = contains(responders_table.classif_data, 'L2/3');
type_oi_neg = contains(non_responders_table.classif_data, 'L2/3');
is_in_exp = ismember(responders_table.experiment, contra);
figure;
scatter(responders_table.L1_cheriff(logical(type_oi_pos)), mm_norm_amps(logical(type_oi_pos)), 65, hex2rgb(responders_table.classif_cols(logical(type_oi_pos))), 'filled');
figure;
scatter(responders_table.L1_cheriff, mm_norm_amps, 65, hex2rgb(responders_table.classif_cols), 'filled');


%% Cell type traces %%

uniqueTypes = unique(responders_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(responders_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);

figure;

for i = 1:size(uniqueTypes)
    if isnan(uniqueTypes{i})
        continue
    end
    ax = subplot(ceil(length(uniqueTypes)/5), 5, i);
    matchType = strcmpi(responders_table.classif_data, uniqueTypes{i});
    matchSlice = ismember(responders_table.experiment', contra)';
    %matchSlice = cell2mat(matchSlice);
    cells = responders_table.snip_data(logical(matchType & matchSlice));
    mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
    mtrace = cell2mat(mtrace)';
    smtrace = -mtrace;
    smtrace = smoothdata(smtrace, 'gaussian', 4);
    
    n = sum(~isnan(mean(smtrace,1)),2);
    upperci = mean(smtrace, 2, 'omitnan')+tinv(0.95,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
    lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.95,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
    
    cols = hex2rgb(uniqueCols{i});

    %plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-1.0, 2.5]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(typeData, 3)))); hold on;
    if size(smtrace,2) > 1
        %plot(upperci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-1, 1]);
        %hold on;
        %plot(lowerci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-1, 1]);
        %hold on;
        %p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.5);
        %p.EdgeColor = cols;
        %hold on;
        %plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(uniqueTypes{i}, ', n=', num2str(size(typeData, 3))));
         % 
        imagesc(smtrace'); clim([-1.0 5]); colormap('pink'); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
        
        %plot(smtrace, 'color', cols, 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 4]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
        %hold on;
        %plot(mean(smtrace, 2, 'omitnan'), 'color', 'black', 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 4]);

        %plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 3.0]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));

    else
        %plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-1, 4]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
    end
end

%% Response properties %%

figure;
matchSlice = find(ismember(response_table.experiment', thalamus)');
rand_idx = randperm(length(matchSlice), 40);
cells = matchSlice(rand_idx);

for i = 1:40
    ax = subplot(8, 5, i);
    cell = cells(i);
    cell = response_table.snip_data(cell);
    mtrace = cellfun(@(x) mean(x, 3), cell, 'UniformOutput', false);
    mtrace = -cell2mat(mtrace)';
    smtrace = smoothdata(mtrace, 'gaussian', 4);

    plot(smtrace, 'color', 'r', 'LineWidth', 1.5); xlim([0 180]); ylim([-1.0, 4]); hold on;
end

matchSlice = find(ismember(responders_table.experiment', thalamus)');
taus = zeros(length(matchSlice), 1);
latency = zeros(length(matchSlice), 1);
rise_times = zeros(length(matchSlice), 1);
for i = 1:length(matchSlice)
    cell = responders_table.snip_data(i);
    mtrace = cellfun(@(x) mean(x, 3), cell, 'UniformOutput', false);
    mtrace = -cell2mat(mtrace)';
    smtrace = smoothdata(mtrace, 'gaussian', 4);
    tau_fit = fit([(1:length(smtrace(43:end)))./400]', smtrace(43:end), 'exp1');
    tau_ci = confint(tau_fit);
    ci_range = (tau_ci(1,2) - tau_ci(2,2))./tau_fit.b;
    if abs(ci_range) > 0.5
        taus(i) = 1000;
        latency(i) = 1000;
        rise_times(i) = 1000;
        continue
    end
    tc = -1/tau_fit.b;
    t_36 = -log(0.3678)*tc;
    taus(i) = t_36;
    d_df = diff(smtrace(39:60));
    for k = 1:(numel(d_df)-2)
        if d_df(k+1) > d_df(k) && d_df(k+2) > d_df(k+1)
            latency(i) = (1+k)/400;
            break
        end
    end
    for j = (k+2):(numel(d_df)-(k+2))
        if d_df(j) < 0
            rise_times(i) = (j-k-2)/400;
            break
        end
    end
    %if i == 30
    %    break
    %end
    %maxValue = max(mtrace(43:63));
    %whereMax = find(mtrace == maxValue);
    %meanValue = mean(mtrace((whereMax-2):(whereMax+2)));
    %poi = find(mtrace(whereMax:120) < meanValue*0.368);
    %if isempty(poi) | length(poi) == 1
    %    taus(i) = -1;
    %else
    %    taus(i) = mean([poi(1) poi(2)])/400;
    %end
end
taus = taus(find(taus < 1000));
rise_times = rise_times(find(rise_times < 1000));
latency = latency(find(latency < 1000));
figure; plot(tau_fit, (1:length(smtrace(43:end)))/400, smtrace(43:end));
figure; boxplot(taus*1000); ylim([0 650]);%xlim([0 150]); % *1000 is in milliseconds

%% Response trace of each cell plotted as a heatmap %%

unique_exps = unique(responders_table.experiment);
allTypes = unique(response_table.classif_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allTypes = allTypes(idx);
responderTypes = categorical(responders_table.classif_data, allTypes);
nonresponderTypes = categorical(non_responders_table.classif_data, allTypes);
types_in_heatmaps = zeros(2,5);

cm = [0 0 0.7; 1 1 1; 0.7 0 0];
cmi = interp1([-1; 0; 1], cm, (-1:0.0001:1));

figure;
for i = 1:numel(unique_exps)
    ax = subplot(2, numel(unique_exps), i);
    %exp_idx = strcmpi(responders_table.experiment, unique_exps{i}); % Choose by experiment
    exp_idx = ismember(responders_table.experiment, contra); % Choose by circuit
    exp_cells = responders_table.snip_data(logical(exp_idx));
    mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
    max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
    max_rows = cell2mat(max_rows);
    norm_rows = cellfun(@(x) x./abs(min(x(PostStim))), mean_rows, 'UniformOutput', false);
    %norm_rows = cellfun(@(x) (x-min(x([1:34 39:end])))/(max(x([1:34 39:end]))-min(x([1:34 39:end]))), mean_rows, 'UniformOutput', false);
    norm_rows = cell2mat(norm_rows);
    % [sort_arr, sort_idx] = sort(max_rows); % Sort just by response
    sorted_types = sort(responderTypes(logical(exp_idx))); % Sort by type and response. Also the next line
    [sort_arr, sort_idx] = sortrows([[grp2idx(responderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
    imagesc(-norm_rows(sort_idx,:)); clim ([-1 1]); colormap(cmi); title(strcat('Responders, n = ', num2str(sum(exp_idx))));
    %yline(find(diff(grp2idx(sorted_types))>=1)+0.5, 'r', 'LineWidth',2); % This finds places where, in the numeric version of the sorted array, the values change. To visually demarcate them
    unique(sorted_types, 'stable')

    ax = subplot(2, numel(unique_exps), i+numel(unique_exps));
    %exp_idx = strcmpi(non_responders_table.experiment, unique_exps{i}); % Choose by experiment
    exp_idx = ismember(non_responders_table.experiment, contra); % Choose by circuit
    exp_cells = non_responders_table.snip_data(logical(exp_idx));
    mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
    max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
    max_rows = cell2mat(max_rows);
    norm_rows = cellfun(@(x) x./abs(max(x([PreStim PostStim PostEPSP]))), mean_rows, 'UniformOutput', false);
    %norm_rows = cellfun(@(x) (x-min(x([1:34 39:end])))/(max(x([1:34 39:end]))-min(x([1:34 39:end]))), mean_rows, 'UniformOutput', false);
    norm_rows = cell2mat(norm_rows);
    %[sort_arr, sort_idx] = sort(max_rows); % Sort just by response
    sorted_types = sort(nonresponderTypes(logical(exp_idx))); % Sort by type and response. Also the next line
    [sort_arr, sort_idx] = sortrows([[grp2idx(nonresponderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
    imagesc(-norm_rows(sort_idx,:)); clim ([-1 1]); colormap(cmi); title(strcat('Non-responders, n = ', num2str(sum(exp_idx))));
    %yline(find(diff(grp2idx(sorted_types))>=1)+0.5, 'r', 'LineWidth',2);
    unique(sorted_types, 'stable')
end
colorbar;

%% Stacked bar of responders vs non responders: %%

allTypes = unique(response_table.classif_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allCols = unique(response_table.classif_cols, 'stable');
allTypes = allTypes(idx);
allCols = allCols(idx);

percents = zeros(numel(allTypes), 3);
stats = zeros(numel(allTypes), 2);

circuit = contra;
for i = 1:numel(allTypes)
    all_in_exp = ismember(response_table.experiment', circuit)';
    resp_in_exp = ismember(responders_table.experiment', circuit)';
    non_resp_in_exp = ismember(non_responders_table.experiment', circuit)';
    matchallType = strcmpi(response_table.classif_data, allTypes{i}) & all_in_exp;
    matchnegType = strcmpi(non_responders_table.classif_data, allTypes{i}) & non_resp_in_exp;
    matchposType = strcmpi(responders_table.classif_data, allTypes{i}) & resp_in_exp;
    all_matches = sum(matchallType);
    neg_matches = sum(matchnegType);
    pos_matches = sum(matchposType);
    all_percent = all_matches/size(response_table(find(all_in_exp),:),1);
    neg_percent = neg_matches/size(non_responders_table(find(non_resp_in_exp),:), 1);
    pos_percent = pos_matches/size(responders_table(find(resp_in_exp), :), 1);
    pval_neg = 2 * min(binocdf(neg_matches, size(non_responders_table(find(non_resp_in_exp),:), 1), all_percent), 1 - binocdf(neg_matches - 1, size(non_responders_table(find(non_resp_in_exp),:), 1), all_percent));
    pval_pos = 2 * min(binocdf(pos_matches, size(responders_table(find(resp_in_exp), :), 1), all_percent), 1 - binocdf(pos_matches - 1, size(responders_table(find(resp_in_exp), :), 1), all_percent));

    percents(i,1) = all_percent;
    percents(i,2) = neg_percent;
    percents(i,3) = pos_percent;

    stats(i,1) = pval_neg;
    stats(i,2) = pval_pos;
end

figure; 
b = bar(percents(:,1:3)', 'stacked', 'FaceColor','flat'); xticklabels({'All responders','Non Responders','Responders'});
ylim([0 1]);
for k = 1:size(percents,1)
    b(k).CData = hex2rgb(allCols{k});
end
stat_summary = {'Under/over represented in:'};
groups = {'Non', 'Resp'};
for j = 1:size(stats,2)
    stat_summary = vertcat(stat_summary, groups{j}, allTypes(find(stats(:,j) < 0.05)));
end
a = annotation('textbox', [0.8 0.5 0.1 0.1], 'String',stat_summary);
a.FontSize = 9;

%% Stacked bar of low and high responders by cell type: %%

allTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allCols = unique(all_table.classif_cols, 'stable');
allTypes = allTypes(idx);
allCols = allCols(idx);

percents = zeros(2, numel(allTypes));
counts = zeros(2, numel(allTypes));
stats = zeros(numel(allTypes),1);
low_resp = all_mm_norm_amps <= 0.02; % norm_amps >= -0.5
high_resp = all_mm_norm_amps > 0.02; %norm_amps < -1.0;

circuit = contra;
for i = 1:numel(allTypes)
    resp_in_exp = ismember(all_table.experiment', circuit)';
    matchlowType = strcmpi(non_responders_table.classif_data, allTypes{i}) & ismember(non_responders_table.experiment, circuit);
    matchhighType = strcmpi(responders_table.classif_data, allTypes{i}) & ismember(responders_table.experiment, circuit);
    low_matches = sum(matchlowType);
    high_matches = sum(matchhighType);
    all_matches = low_matches + high_matches;
    low_percent = low_matches/all_matches;
    high_percent = high_matches/all_matches;
    pval_high_to_low = binocdf(high_matches, all_matches, high_percent);

    percents(1,i) = low_percent;
    percents(2,i) = high_percent;
    counts(1,i) = low_matches;
    counts(2,i) = high_matches;

    stats(i) = pval_high_to_low;
end

figure('Position', [10 10 1000 500]);
h1 = axes;
b1 = bar(percents(2,:), 'FaceColor', 'flat');
b1(1).CData = hex2rgb(allCols);
set(h1,'Ydir','normal');
xticks(1:numel(allTypes)); xticklabels(allTypes); xtickangle(-45);

figure('Position', [10 10 1000 500]);
h1 = axes;
b1 = bar(counts(1,:), 'FaceColor', 'flat');
b1(1).CData = hex2rgb(allCols);
set(h1,'Ydir','normal'); hold on;
h2 = gca;
b2 = bar(-counts(2,:), 'FaceColor', 'flat'); 
set(h2, 'Ydir', 'reverse');
b2(1).CData = hex2rgb(allCols);
xticks(1:numel(allTypes)); xticklabels(allTypes); xtickangle(-45);
%ylim([-30 30]);

%% Stacked bar cell types by slice (to check for appropriate normalization): %%

allTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allCols = unique(all_table.classif_cols, 'stable');
allTypes = allTypes(idx);
allCols = allCols(idx);
allSlices = unique(responders_table.slice_data);
circuits = cell(numel(all_mm_norm_amps),1);
table_to_use = all_table(1:numel(all_mm_norm_amps), :);
for i = 1:numel(circuits)
    if ismember(table_to_use.experiment(i), thalamus)
        circuits(i) = {'thalamus'};
    else
        circuits(i) = {'contra'};
    end
end
allCircuits = unique(circuits);
percents = zeros(numel(allCircuits), numel(allTypes), 1);

circuit = thalamus;
for i = 1:numel(allCircuits)
    %who = logical(ismember(responders_table.experiment', circuit)' & responders_table.slice_data == i);
    who = logical(ismember(circuits, allCircuits{i})');
    resp_in_exp = all_table(who,:);
    for j = 1:numel(allTypes)
            matchposType = strcmpi(resp_in_exp.classif_data, allTypes{j});
            pos_matches = sum(matchposType);
            pos_percent = pos_matches/size(resp_in_exp, 1);
        
            percents(i,j) = pos_percent;
    end
end
    
figure; 
b = bar(percents, 'stacked', 'FaceColor','flat');
ylim([0 1]);
for k = 1:size(percents,2)
    b(k).CData = hex2rgb(allCols{k});
end

%% Stacked bar of responders vs non responders by laminar "bin": %%

edges = [-1, -0.9, -0.7, -0.55, -0.4, -0.25, -0.075, 0]; % Laminar bins based on subjective boundaries. First bin is L6b, last bin is L1
binned_depth_all = discretize(all_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'LL5b', 'UL5b', 'L4/5a', 'L2/3', 'L1'});
binned_depth_pos = discretize(responders_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'LL5b', 'UL5b', 'L4/5a', 'L2/3', 'L1'});

allTypes = unique(response_table.classif_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allCols = unique(response_table.classif_cols, 'stable');
allTypes = allTypes(idx);
allCols = allCols(idx);

percents = zeros(numel(edges)-1, numel(allTypes), 2);
stats = zeros(numel(edges)-1, numel(allTypes));

for j = 1:numel(edges)-1
    all_in_bin = response_table(find(binned_depth_all == j),:);
    pos_in_bin = responders_table(find(binned_depth_pos == j & norm_amps >= -0.5),:); %& norm_amps < -0.5
    if size(pos_in_bin, 1) == 0
        if size(all_in_bin) == 0
            percents(j,i,1) = 0;
            percents(j,i,2) = 0;
            
            stats(j,i) = 2;
        else
            for i = 1:numel(allTypes)
                matchallType = strcmpi(all_in_bin.classif_data, allTypes{i});
                all_matches = sum(matchallType);
                all_percent = all_matches/size(all_in_bin, 1);
            
                percents(j,i,1) = all_percent;
                percents(j,i,2) = 0;
            
                stats(j,i) = 2;
            end
        end
    else
        for i = 1:numel(allTypes)
            matchallType = strcmpi(all_in_bin.classif_data, allTypes{i});
            matchposType = strcmpi(pos_in_bin.classif_data, allTypes{i});
            all_matches = sum(matchallType);
            pos_matches = sum(matchposType);
            all_percent = all_matches/size(all_in_bin, 1);
            pos_percent = pos_matches/size(pos_in_bin, 1);
            pval_pos = 2 * min(binocdf(pos_matches, size(pos_in_bin, 1), all_percent), 1 - binocdf(pos_matches - 1, size(pos_in_bin, 1), all_percent));
        
            percents(j,i,1) = all_percent;
            percents(j,i,2) = pos_percent;
        
            stats(j,i) = pval_pos;
        end
    end
end

figure;
for i = 1:size(percents,3)
   ax = subplot(size(percents,3), 1, i);
   b = bar(percents(:,:,i), 'stacked', 'FaceColor','flat'); hold on;
   ylim([0 1]);
   for k = 1:size(percents,2)
       b(k).CData = hex2rgb(allCols{k});
   end
end
stat_summary = {'Under/over represented in responding:'};
for j = 1:size(stats,1)
    stat_summary = vertcat(stat_summary, strcat('Bin ', num2str(j)), allTypes(find(stats(j,:) < 0.05)));
end
a = annotation('textbox', [0.8 0.3 0.1 0.1], 'String',stat_summary);
a.FontSize = 9;

%% Violin and bar plots for comparing cell type responses across circuits

pos_mm_norm_amps = all_mm_norm_amps(find(all_mm_norm_amps > -1));
log_mm_norm_amps = log2(pos_mm_norm_amps+1);
%log_mm_norm_amps = max(-0.5, log_mm_norm_amps);

contra_grouping = categorical(all_table.subclass_data(logical(ismember(all_table.experiment, contra) & all_mm_norm_amps > -1)), ...
    allTypes);
thalamus_grouping = categorical(all_table.subclass_data(logical(ismember(all_table.experiment, thalamus) & all_mm_norm_amps > -1)), ...
    allTypes);
[c1,cf1] = kde(log_mm_norm_amps(logical(ismember(all_table.experiment, contra))), Bandwidth=0.05);
figure;
violinplot(contra_grouping, log_mm_norm_amps(logical(ismember(all_table.experiment, contra))), ...
    DensityDirection='negative', DensityScale='width');
yline(1); hold on; %ylim([-1 12]);
violinplot(thalamus_grouping, log_mm_norm_amps(logical(ismember(all_table.experiment, thalamus))), ...
    DensityDirection='positive', DensityScale='width');
yline(1);

grouping = categorical(all_table.subclass_data(find(all_mm_norm_amps > -1)), allTypes);
color_grouping = ismember(all_table.experiment(find(all_mm_norm_amps > -1)), thalamus);
figure;
boxchart(grouping, pos_mm_norm_amps, "GroupByColor", color_grouping, 'Notch', 'on'); hold on;
swarmchart(grouping(~logical(color_grouping)), pos_mm_norm_amps(~logical(color_grouping)), 'b.', 'XJitterWidth', 1); hold on;
swarmchart(grouping(logical(color_grouping)), pos_mm_norm_amps(logical(color_grouping)), 'r.', 'XJitterWidth', 0.25);
yline(1); ylabel('Median-normalized response amplitude (au)');

figure;
swarmchart(grouping, log(abs(all_mm_norm_amps)), '.'); ylim([-5 50]);

%% Violin/boxplots for responses by depth in cortex

edges = [-1, -0.9, -0.7, -0.55, -0.4, -0.25, -0.075, 0]; % Laminar bins based on subjective boundaries. First bin is L6b, last bin is L1
binned_depth_all = discretize(all_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'LL5b', 'UL5b', 'L4/5a', 'L2/3', 'L1'});
binned_depth_pos = discretize(responders_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'LL5b', 'UL5b', 'L4/5a', 'L2/3', 'L1'});
type_grouping = categorical(all_table.subclass_data, allTypes);

% Trying this by hand:
grouping_to_plot = binned_depth_all; %type_grouping
swarm_colors = repmat([0 0 0],numel(all_mm_norm_amps),1);
swarm_colors(1:numel(mm_norm_amps),:) = repmat([1 0 0],numel(mm_norm_amps),1);
swarm_colors((numel(mm_norm_amps)+1):end,:) = repmat([0.5 0.5 0.5],numel(all_mm_norm_amps)-numel(mm_norm_amps),1);
exper = thalamus;
include = logical(ismember(all_table.experiment, exper));
resp_or_not = zeros(numel(all_mm_norm_amps),1);
resp_or_not(1:numel(mm_norm_amps)) = ones(numel(mm_norm_amps),1);
figure;
violinplot(grouping_to_plot(include), log2(all_mm_norm_amps(include)+1), DensityScale='width', DensityDirection='negative'); 
ylim([-1.5 5]); yline(1); hold on;
%swarmchart(grouping_to_plot(include), log2(all_mm_norm_amps(include)+1), [],...
%    swarm_colors(find(include),:), 'XJitterWidth', 0.25); hold on;
count = 1;
for i = categories(grouping_to_plot)'
    in_bin_resp = log2(all_mm_norm_amps(include & grouping_to_plot == i & resp_or_not)+1);
    in_bin_non = log2(all_mm_norm_amps(include & grouping_to_plot == i & ~resp_or_not)+1);
    resp_med = median(in_bin_resp);
    non_med = median(in_bin_non);
    plot(count, resp_med, 'k', 'Marker', 'diamond', 'MarkerSize', 12, 'MarkerFaceColor', 'k'); hold on;
    plot(count, non_med, 'b', 'Marker', 'diamond', 'MarkerSize', 12, 'MarkerFaceColor', 'b'); hold on;
    count = count+1;
end

% % Using grpandplot function
% circuits = cell(numel(all_mm_norm_amps),1);
% table_to_use = all_table(1:numel(all_mm_norm_amps), :);
% for i = 1:numel(circuits)
%     if ismember(table_to_use.experiment(i), thalamus)
%         circuits(i) = {'thalamus'};
%     else
%         circuits(i) = {'contra'};
%     end
% end
% 
% % All cells
% parent_fig = figure;
% grpandplot(table(all_mm_norm_amps, binned_depth_all, circuits, ...
%     'VariableNames', {'data', 'depth', 'circuit'}), ...
%     'data', xFactor='depth', cFactor='circuit', showVln = true, showBox=false, ...
%     pntSize = 5, parent=parent_fig); ylim([-1.5 25]); hold on;
% yline(1); ylabel('Median-normalized response amplitude (au)');
% 
% % Non-responders
% grpandplot(table(all_mm_norm_amps((numel(mm_norm_amps)+1):end), binned_depth_all((numel(mm_norm_amps)+1):end), circuits((numel(mm_norm_amps)+1):end), ...
%     'VariableNames', {'data', 'depth', 'circuit'}), ...
%     'data', xFactor='depth', cFactor='circuit', showVln = false, showBox=false, ...
%     pntSize = 5, pntFillC=[0 0 0], parent=parent_fig); ylim([-1.5 25]);
% yline(1); ylabel('Median-normalized response amplitude (au)'); %ylim([-5 100]);
% 
% % Responders
% grpandplot(table(all_mm_norm_amps(1:numel(mm_norm_amps)), binned_depth_all(1:numel(mm_norm_amps)), circuits(1:numel(mm_norm_amps)), ...
%     'VariableNames', {'data', 'depth', 'circuit'}), ...
%     'data', xFactor='depth', cFactor='circuit', showVln = true, showBox=false, ...
%     pntSize = 5); ylim([-1.5 25]);
% yline(1); ylabel('Median-normalized response amplitude (au)'); %ylim([-5 100]);

%% Response amplitudes relative to columnar cheriff %%

allTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allCols = unique(responders_table.classif_cols, 'stable');
allTypes = allTypes(idx);
grouping = categorical(all_table.classif_data, allTypes);
color_grouping = ismember(all_table.experiment, contra);
combined_color_grouping = cat(1,grouping(color_grouping), grouping(color_grouping));
combined_amps = cat(1, all_mm_norm_amps(color_grouping), all_mm_norm_amps(color_grouping)./all_table.cheriff_data(color_grouping));
cheriff_group = cat(2, repelem(1, numel(combined_amps)/2), repelem(2, numel(combined_amps)/2));

figure('Position', [10 10 1000 900]);
ax = subplot(2, 1, 1);
boxchart(combined_color_grouping, combined_amps, 'GroupByColor', cheriff_group);
title('Responses'); ylabel('Median-normalized response amp');
%ylim([-1 12]);
ax = subplot(2, 1, 2);
boxchart(grouping(color_grouping), all_mm_norm_amps(color_grouping)./all_table.cheriff_data(color_grouping));
title('Responses normalized to columnar CheRiff signal'); ylabel('Median-normalized response amp');
%ylim([-1 12]);

figure;
scatter(responders_table.cheriff_data, mm_norm_amps, 65, hex2rgb(responders_table.classif_cols), 'filled');
ylim([0 10]);
xlabel('Columnar Cheriff value, slice normalized');
ylabel('Response amplitude');


%% Response amplitudes plotted by slice depth (axial to imaging) %%

figure;
res_plot = scatter(responders_table.slice_depth, mm_norm_amps, 65, hex2rgb(responders_table.classif_cols), 'filled');
xlabel('Axial depth in slice');
ylabel('Response amplitude');

%% Cell type traces, with 90% confidence intervals:

uniqueTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(all_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);

boot_ci_means = NaN(numel(uniqueTypes),1);
boot_ci_sem = NaN(numel(uniqueTypes),1);
boot_ci_n = NaN(numel(uniqueTypes), 1);
whole_boot_ci = NaN(size(uniqueTypes,1),50000);

figure('Position', [10 10 1000 900]);

for i = 1:size(uniqueTypes)
    if isnan(uniqueTypes{i})
        continue
    end
    ax = subplot(ceil(length(uniqueTypes)/5), 5, i);
    matchType = strcmpi(all_table.classif_data, uniqueTypes{i});
    matchSlice = ismember(all_table.experiment', thalamus)';
    %matchSlice = strcmpi(all_table.experiment, '0730');
    cells = all_table.snip_data(logical(matchType & matchSlice));
    num_cells = sum(logical(matchType & matchSlice));
    % matchType_resp = strcmpi(responders_table.subclass_data, uniqueTypes{i});
    % matchSlice_resp = ismember(responders_table.experiment', contra)';
    % num_cells_resp = sum(logical(matchType_resp & matchSlice_resp));
    % disp(strcat(uniqueTypes{i},' _resp_prop = ',num2str(num_cells_resp/num_cells)));
    mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
    mtrace = cell2mat(mtrace)';
    mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
    smtrace = -mtrace;
    smtrace = smoothdata(smtrace, 'gaussian', 8);
    if size(smtrace,2) == 0
        continue
    end
    
    cols = hex2rgb(uniqueCols{i});
    
    if size(smtrace,2) > 1
        n = sum(~isnan(mean(smtrace,1)),2);
        [ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
        upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
        lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
        %upperci = smoothdata(ci(2,:)', 'gaussian', 8);
        %lowerci = smoothdata(ci(1,:)', 'gaussian', 8);
        %upperci = mean(smtrace, 2, 'omitnan')+tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
        %lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
        %plot(upperci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]);
        %hold on;
        %plot(lowerci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-0.5, 6]);
        %hold on;
        p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.3); xlim([0 180]);
        p.EdgeColor = cols;
        hold on;
        plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2)))); ylim([-0.5, 2]);
        disp(strcat(uniqueTypes{i},' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
        %max_idx = find(mean(smtrace, 2, 'omitnan') == max(mean(smtrace, 2, 'omitnan')));
        disp(strcat('se = ', num2str(std(bootstat(:,PostStim),0,'all'))));
        boot_ci_means(i) = mean(smtrace(PostStim,:), 'all', 'omitnan');
        boot_ci_sem(i) = std(bootstat(:,PostStim),0,'all');
        boot_ci_n(i) = size(smtrace, 2);
        yline(0, '--');
    else
        plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-0.5, 2]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
        disp(strcat(uniqueTypes{i},' = ',num2str(max(smtrace))));
        yline(0, '--');
    end
end

%% Excit or Inhib cell type traces, with 90% confidence intervals, by depth:

uniqueTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(all_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);
edges = [-1, -0.9, -0.7, -0.4, -0.25, -0.075, 0]; % Laminar bins based on subjective boundaries. First bin is L6b, last bin is L1
binned_depth_all = discretize(all_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'L5b', 'L4/5a', 'L2/3', 'L1'});
all_depths = flip(unique(binned_depth_all));
inhib_types = uniqueTypes(contains(uniqueTypes, 'Gaba'));
excit_types = uniqueTypes(contains(uniqueTypes, 'Glut')); %& contains(uniqueTypes, 'L5 IT')

types_oi = inhib_types;
figure('Position', [10 10 1000 900]);
count = 1;
for i = 1:numel(types_oi)
    for k = 1:(numel(all_depths))
        ax1 = subplot(numel(types_oi), numel(all_depths), count);
        matchType = strcmpi(all_table.classif_data, types_oi{i});
        matchSlice = ismember(all_table.experiment', contra)';
        matchDepth = binned_depth_all == all_depths(k);
        cells = all_table.snip_data(logical(matchType & matchSlice & matchDepth));
        if ~numel(cells)
            count = count + 1;
            continue
        else
            mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
            mtrace = cell2mat(mtrace)';
            mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
            smtrace = -mtrace;
            smtrace = smoothdata(smtrace, 'gaussian', 8);
            if size(smtrace,2) == 0
                continue
            end
        
            cols = hex2rgb(uniqueCols{i});
            
            if size(smtrace,2) > 1
                n = sum(~isnan(mean(smtrace,1)),2);
                [ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
                upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
                lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
                %upperci = smoothdata(ci(2,:)', 'gaussian', 8);
                %lowerci = smoothdata(ci(1,:)', 'gaussian', 8);
                %upperci = mean(smtrace, 2, 'omitnan')+tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
                %lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
                plot(upperci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]);
                hold on;
                plot(lowerci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-0.5, 7]);
                hold on;
                p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.5);
                p.EdgeColor = cols;
                hold on;
                plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(types_oi{i}, ', n=', num2str(size(smtrace, 2)), {sprintf('\n')}, string(all_depths(k))));

                %plot(smtrace, 'color', cols, 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 5]); title(strcat(excit_types{i}, ', n=', num2str(size(smtrace, 2)), {sprintf('\n')}, string(all_depths(k))));
                %hold on;
                %plot(mean(smtrace, 2, 'omitnan'), 'color', 'black', 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 4]);
                disp(strcat(types_oi{i}, string(all_depths(k)),' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
                %max_idx = find(mean(smtrace, 2, 'omitnan') == max(mean(smtrace, 2, 'omitnan')));
                disp(strcat('se = ', num2str(std(bootstat(:,PostStim),0,'all'))));
                yline(0, '--');
            else
                plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-0.5, 7]); title(strcat(types_oi{i}, ', n=', num2str(size(smtrace, 2)), {sprintf('\n')}, string(all_depths(k))));
                disp(strcat(types_oi{i}, string(all_depths(k))));
            end
            count = count + 1;
        end
    end
end

%% Response traces by depth
edges = [-1, -0.9, -0.7, -0.4, -0.25, -0.075, 0]; % Laminar bins based on subjective boundaries. First bin is L6b, last bin is L1
binned_depth_all = discretize(all_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'L5b', 'L4/5a', 'L2/3', 'L1'});
binned_depth_pos = discretize(responders_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'L5b', 'L4/5a', 'L2/3', 'L1'});
all_depths = flip(unique(binned_depth_all));

figure('Position', [10 10 1000 200]);
count = 1;
for i = 1:(numel(all_depths))
    ax1 = subplot(1, numel(all_depths), count);
    %matchSlice = ismember(all_table.experiment', thalamus)';
    matchSlice = strcmpi(all_table.experiment, '0814');
    matchDepth = binned_depth_all == all_depths(i);
    cells = all_table.snip_data(logical(matchSlice & matchDepth));
    if ~numel(cells)
        count = count + 1;
        continue
    else
        mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
        mtrace = cell2mat(mtrace)';
        mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
        smtrace = -mtrace;
        smtrace = smoothdata(smtrace, 'gaussian', 8);
        if size(smtrace,2) == 0
            continue
        end
    
        cols = 'b';
        
        if size(smtrace,2) > 1
            n = sum(~isnan(mean(smtrace,1)),2);
            [ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
            upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
            lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
            %upperci = smoothdata(ci(2,:)', 'gaussian', 8);
            %lowerci = smoothdata(ci(1,:)', 'gaussian', 8);
            %upperci = mean(smtrace, 2, 'omitnan')+tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
            %lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
            plot(upperci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]);
            hold on;
            plot(lowerci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-0.5, 2.1]);
            hold on;
            p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.5);
            p.EdgeColor = cols;
            hold on;
            plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(string(all_depths(i)), ', n=', num2str(size(smtrace, 2))));
            disp(strcat(uniqueTypes{i},' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
            %max_idx = find(mean(smtrace, 2, 'omitnan') == max(mean(smtrace, 2, 'omitnan')));
            disp(strcat('se = ', num2str(std(bootstat(:,PostStim),0,'all'))));
            %plot(smtrace, 'color', cols, 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 5]); title(strcat(excit_types{i}, ', n=', num2str(size(smtrace, 2)), {sprintf('\n')}, string(all_depths(k))));
            %hold on;
            %plot(mean(smtrace, 2, 'omitnan'), 'color', 'black', 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 4]);

        else
            plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-0.5, 2.1]); title(strcat(string(all_depths(i)), ', n=', num2str(size(smtrace, 2))));
        end
        count = count + 1;
    end
end

%% Statistical tests across cell type traces

uniqueTypes(7)

one = 6;
two = 11;

mean_one = contra_ci_means(16);
mean_two = thalamus_ci_means(16);
sem_one = contra_ci_sem(16);
sem_two = thalamus_ci_sem(16);
n_one = contra_ci_n(16);
n_two = thalamus_ci_n(16);

% t stat for comparing two groups:
sp = sqrt(((n_one - 1)*sem_one^2 + (n_two - 1)*sem_two^2)/(n_one + n_two - 2));
t = (mean_one-mean_two)/(sp * sqrt(1/n_one + 1/n_two));
t
(n_one + n_two - 2)

% t stat for comparing whole circuits:
total_n = size(all_table(logical(ismember(all_table.experiment', contra)),:),1);
total_amps = cellfun(@(x) mean(x(:,PostStim,:), [2 3]), all_table.snip_data(logical(ismember(all_table.experiment', contra))), 'UniformOutput', false);
total_amps = cell2mat(total_amps)';
total_means = mean(total_amps);
total_std = std(total_amps, 0, 'all');

% ANOVA for exp vs cell type
g1 = vertcat(uniqueTypes,uniqueTypes);
g2 = vertcat(repmat('0730', numel(uniqueTypes),1),repmat('0802',numel(uniqueTypes),1));
stat_struct = vertcat(thalamus_ci_means(:,1), thalamus_ci_means(:,2));
[~,~,stats] = anovan(stat_struct, {g1,g2});
[results,~,~,gnames] = multcompare(stats, 'Dimension', [1 2]);

% Standard error of the bootstrapped ratios of means between normalized cell types
ratios = (whole_boot_ci(two,:))./(whole_boot_ci(one,:));
std(ratios, 0, 2, 'omitnan')

%% Response amplitude scatter to compare cell type responses across circuits

uniqueTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(all_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);

figure;
for i = 1:(numel(uniqueTypes))

    matchSlice = ismember(all_table.experiment', contra)';
    matchType = strcmpi(all_table.classif_data, uniqueTypes{i});
    contra_cells = all_mm_norm_amps(logical(matchSlice & matchType));
    mean_contra_cells = mean(contra_cells);
    sem_contra_cells = std(contra_cells)/sqrt(numel(contra_cells));
    
    matchSlice = ismember(all_table.experiment', thalamus)';
    thalamus_cells = all_mm_norm_amps(logical(matchSlice & matchType));
    mean_thalamus_cells = mean(thalamus_cells);
    sem_thalamus_cells = std(thalamus_cells)/sqrt(numel(thalamus_cells));
    
    cols = hex2rgb(uniqueCols{i});
    
    scatter(mean_contra_cells, mean_thalamus_cells, 50, cols, 'filled'); hold on;
    ex = errorbar(mean_contra_cells, mean_thalamus_cells, sem_contra_cells, 'horizontal'); hold on;
    ey = errorbar(mean_contra_cells, mean_thalamus_cells, sem_thalamus_cells, 'vertical');
    set(ex, 'color', cols, 'LineWidth', 1);
    set(ey, 'color', cols, 'LineWidth', 1);
    text(mean_contra_cells, mean_thalamus_cells, strcat(uniqueTypes{i}, ', n=', num2str(numel(contra_cells)), '; n=', num2str(numel(thalamus_cells))));
    
end
xlim([-0.5 4.5]); ylim([-0.5 4.5]);
line(xlim, ylim);

%% Response heatmaps for Mike:

allTypes = unique(response_table.subclass_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allTypes = allTypes(idx);
responderTypes = categorical(responders_table.subclass_data, allTypes);
nonresponderTypes = categorical(non_responders_table.subclass_data, allTypes);
types_in_heatmaps = zeros(2,5);

cm = [0 0 0.7; 1 1 1; 0.7 0 0];
cmi = interp1([-1; 0; 1], cm, (-1:0.0001:1));

figure;

ax = subplot(1, 2, 1);
exp_idx = ismember(responders_table.experiment, thalamus);
exp_cells = responders_table.snip_data(exp_idx);
mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
max_rows = cell2mat(max_rows);
norm_rows = cellfun(@(x) x./abs(min(x(PostStim))), mean_rows, 'UniformOutput', false);
%norm_rows = cellfun(@(x) (x-min(x([1:34 39:end])))/(max(x([1:34 39:end]))-min(x([1:34 39:end]))), mean_rows, 'UniformOutput', false);
norm_rows = cell2mat(norm_rows);
sorted_types = sort(responderTypes(logical(exp_idx)));
[sort_arr, sort_idx] = sortrows([[grp2idx(responderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
imagesc(-norm_rows(randperm(size(norm_rows,1)),:));
yline(find(diff(grp2idx(sorted_types))>=1)+0.5, 'w', 'LineWidth',2); % This finds places where, in the numeric version of the sorted array, the values change. To visually demarcate them
clim ([-1 1]); 
colormap(cmi); 
title(strcat('Responders, n = ', num2str(size(exp_cells,1))));
colorbar;

ax = subplot(1, 2, 2);
exp_idx = ismember(non_responders_table.experiment, thalamus);
exp_cells = non_responders_table.snip_data(exp_idx);
mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
max_rows = cell2mat(max_rows);
norm_rows = cellfun(@(x) x./abs(max(x(PostStim))), mean_rows, 'UniformOutput', false);
%norm_rows = cellfun(@(x) (x-min(x([1:34 39:end])))/(max(x([1:34 39:end]))-min(x([1:34 39:end]))), mean_rows, 'UniformOutput', false);
norm_rows = cell2mat(norm_rows);
sorted_types = sort(nonresponderTypes(logical(exp_idx)));
[sort_arr, sort_idx] = sortrows([[grp2idx(nonresponderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
imagesc(-norm_rows(randperm(size(norm_rows,1)),:));
yline(find(diff(grp2idx(sorted_types))>=1)+0.5, 'w', 'LineWidth',2); % This finds places where, in the numeric version of the sorted array, the values change. To visually demarcate them
clim ([-1 1]); 
colormap(cmi); 
title(strcat('Non-responders, n = ', num2str(size(exp_cells,1))));
colorbar;


%% Giant array of cell type responses

matchExp = ismember(all_table.experiment', {'0811'})';
matchSlice = ismember(all_table.slice_data, [1 2 3 4 5]);
cells = all_table.snip_data(logical(matchExp & matchSlice));
respond_bool = is_responding(logical(matchExp & matchSlice));
corr_amps = cellfun(@(x) mean(x(1,PostStim,:), 'all')-mean(x(1,PreStim,:), 'all'), cells, 'UniformOutput', false);
amp_bins = discretize(-cell2mat(corr_amps), 5);

cmap = [
    0.5, 0.5, 0.5;  % Grey (Bin 1)
    0.8, 0.5, 0.5;  % Light red (pink) (Bin 2)
    0.9, 0.4, 0.4;  % Soft red (Bin 3)
    0.9, 0.2, 0.2;  % Deeper red (Bin 4)
    0.8, 0, 0;      % Saturated red (Bin 5)
];

idx = randsample(size(cells,1), 209, false);

figure('Position', [10 10 1000 1000]);

for i = 1:numel(idx)
    idx_use = idx(i);

    ax = subplot(11, 19, i);
    cell = cells(idx_use);
    mtrace = mean(cell2mat(cell),3);
    smtrace = -mtrace;
    smtrace(Stim) = zeros(numel(Stim),1);
    smtrace = smoothdata(smtrace, 'gaussian', 6);
    cols = cmap(amp_bins(idx_use),:);

    % if respond_bool(idx_use)
    %     cols = cmap(amp_bins(idx_use),:);
    % else
    %     cols = [0.5 0.5 0.5];
    % end

    plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-1.0, 5]);
    yline(0, '--');
    hold on;
    
end
colormap(cmap);
colorbar;

%% Response trace of each cell in a particular dataset plotted as a heatmap %%

unique_exps = flip(unique(responders_table.experiment));
allTypes = unique(response_table.subclass_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allTypes = allTypes(idx);
responderTypes = categorical(responders_table.subclass_data, allTypes);
nonresponderTypes = categorical(non_responders_table.subclass_data, allTypes);
types_in_heatmaps = zeros(2,5);

cm = [0 0 0.7; 1 1 1; 0.7 0 0];
cmi = interp1([-1; 0; 1], cm, (-1:0.0001:1));

figure;
for i = 1:numel(unique_exps)
    ax = subplot(2, numel(unique_exps), i);
    exp_idx = strcmpi(responders_table.experiment, unique_exps{i});
    %exp_idx = ismember(responders_table.experiment, lo);
    exp_cells = responders_table.snip_data(logical(exp_idx));
    mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
    max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
    max_rows = cell2mat(max_rows);
    norm_rows = cellfun(@(x) x./abs(min(x(PostStim))), mean_rows, 'UniformOutput', false);
    %norm_rows = cellfun(@(x) (x-min(x([1:34 39:end])))/(max(x([1:34 39:end]))-min(x([1:34 39:end]))), mean_rows, 'UniformOutput', false);
    norm_rows = cell2mat(norm_rows);
    [sort_arr, sort_idx] = sort(max_rows);
    %sorted_types = sort(responderTypes(logical(exp_idx)));
    %[sort_arr, sort_idx] = sortrows([[grp2idx(responderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
    imagesc(-norm_rows(sort_idx,:)); clim ([-1 1]); colormap(cmi); title(strcat('Responders, n = ', num2str(sum(exp_idx))));
    %yline(find(diff(grp2idx(sorted_types))>=1)+0.5, 'r', 'LineWidth',2); % This finds places where, in the numeric version of the sorted array, the values change. To visually demarcate them
    %unique(sorted_types, 'stable')

    ax = subplot(2, numel(unique_exps), i+numel(unique_exps));
    exp_idx = strcmpi(non_responders_table.experiment, unique_exps{i});
    %exp_idx = ismember(non_responders_table.experiment, lo);
    exp_cells = non_responders_table.snip_data(logical(exp_idx));
    mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
    max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
    max_rows = cell2mat(max_rows);
    norm_rows = cellfun(@(x) x./abs(max(x(PostStim))), mean_rows, 'UniformOutput', false);
    %norm_rows = cellfun(@(x) (x-min(x([1:34 39:end])))/(max(x([1:34 39:end]))-min(x([1:34 39:end]))), mean_rows, 'UniformOutput', false);
    norm_rows = cell2mat(norm_rows);
    [sort_arr, sort_idx] = sort(max_rows);
    %sorted_types = sort(nonresponderTypes(logical(exp_idx)));
    %[sort_arr, sort_idx] = sortrows([[grp2idx(nonresponderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
    imagesc(-norm_rows(sort_idx,:)); clim ([-1 1]); colormap(cmi); title(strcat('Non-responders, n = ', num2str(sum(exp_idx))));
    %yline(find(diff(grp2idx(sorted_types))>=1)+0.5, 'r', 'LineWidth',2);
    %unique(sorted_types, 'stable')
end
colorbar;

%% Average response trace of binned responders in lo vs hi TTX

lo_cells = all_table(strcmp(all_table.experiment, '0811'),:);
lo_resps = all_mm_norm_amps(strcmp(all_table.experiment, '0811'));
hi_cells = all_table(strcmp(all_table.experiment, '0730'),:);
hi_resps = all_mm_norm_amps(strcmp(all_table.experiment, '0730'));
N_bins = 6;
lo_q = quantile(lo_resps, 0:1/N_bins:1);
hi_q = quantile(hi_resps, 0:1/N_bins:1);
[disc_lo, N_lo] = discretize(lo_resps, lo_q);
[disc_hi, N_hi] = discretize(hi_resps, hi_q);

figure;
for i = 1:N_bins
    ax = subplot(2, N_bins, i);
    lo_to_use = lo_cells.snip_data(disc_lo == i);
    lo_mtrace = cellfun(@(x) mean(x, 3), lo_to_use, 'UniformOutput', false);
    lo_mtrace = -cell2mat(lo_mtrace)';
    lo_mtrace = smoothdata(lo_mtrace, 1, 'gaussian', 4);
    lo_mtrace(Stim,:) = zeros(numel(Stim),size(lo_mtrace,2));
    lo_smtrace = mean(lo_mtrace,2);
    lo_smtrace(Stim) = zeros(numel(Stim),1);
    ci = bootci(500, {@(x) mean(x,1,'omitnan'), lo_mtrace'}, Options=statset(UseParallel=true));
    upperci = ci(2,:)';
    sm_upperci = smoothdata(upperci, 'gaussian', 4);
    lowerci = ci(1,:)';
    sm_lowerci = smoothdata(lowerci, 'gaussian', 4);
    %upperci = mean(smtrace, 2, 'omitnan')+tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
    %lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
    plot(sm_upperci, 'LineWidth', 0.1);
    hold on;
    plot(sm_lowerci,'LineWidth', 0.1);
    hold on;
    p = fill([(1:size(sm_lowerci,1)).'; flip(1:size(sm_upperci,1)).'],[sm_lowerci; flip(sm_upperci)], 'b', 'FaceAlpha', 0.5);
    p.EdgeColor = 'b';
    hold on;

    plot(lo_smtrace); yline(0); ylim([-1 12]); xlim([0 200]); hold on;
    
    ax = subplot(2, N_bins, i+N_bins);
    hi_to_use = hi_cells.snip_data(disc_hi == i);
    hi_mtrace = cellfun(@(x) mean(x, 3), hi_to_use, 'UniformOutput', false);
    hi_mtrace = -cell2mat(hi_mtrace)';
    hi_mtrace = smoothdata(hi_mtrace, 1, 'gaussian', 4);
    hi_mtrace(Stim,:) = zeros(numel(Stim),size(hi_mtrace,2));
    hi_smtrace = mean(hi_mtrace,2);
    ci = bootci(500, {@(x) mean(x,1,'omitnan'), hi_mtrace'}, Options=statset(UseParallel=true));
    upperci = ci(2,:)';
    sm_upperci = smoothdata(upperci, 'gaussian', 4);
    lowerci = ci(1,:)';
    sm_lowerci = smoothdata(lowerci, 'gaussian', 4);
    %upperci = mean(smtrace, 2, 'omitnan')+tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
    %lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
    plot(sm_upperci, 'LineWidth', 0.1);
    hold on;
    plot(sm_lowerci,'LineWidth', 0.1);
    hold on;
    p = fill([(1:size(sm_lowerci,1)).'; flip(1:size(sm_upperci,1)).'],[sm_lowerci; flip(sm_upperci)], 'b', 'FaceAlpha', 0.5);
    p.EdgeColor = 'b';
    hold on;

    plot(hi_smtrace); yline(0); ylim([-1 12]); xlim([0 200]); hold on;
end

%% Circuit PSP sizes:


matchSlice = ismember(all_table.experiment', contra)';
matchClass = contains(all_table.subclass_data, 'Glut');
matchType = contains(all_table.subclass_data, 'IT');
matchDepth = all_table.depth_data < -0.5;
cells = all_table.snip_data(logical(matchSlice & matchClass & ~matchDepth));
mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
mtrace = cell2mat(mtrace)';
mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
smtrace = -mtrace;
smtrace = smoothdata(smtrace, 'gaussian', 8);

n = sum(~isnan(mean(smtrace,1)),2);
[ci, bootstat] = bootci(1000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);

figure;
p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.3); xlim([0 180]);
p.EdgeColor = cols;
hold on;
plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat('Deep cells, n=', num2str(size(smtrace, 2)))); ylim([-0.5, 2]);
yline(0, '--');

disp(strcat('thalamus',' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
disp(strcat('se = ', num2str(std(bootstat(:,PostStim),0,'all'))));

%matchSlice = ismember(all_table.experiment', contra)';
matchSlice = strcmpi(all_table.experiment, '0814');
cells = all_table.snip_data(logical(matchSlice));
mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
mtrace = cell2mat(mtrace)';
mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
smtrace = -mtrace;
smtrace = smoothdata(smtrace, 'gaussian', 8);

n = sum(~isnan(mean(smtrace,1)),2);
[ci, bootstat] = bootci(1000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);

disp(strcat('contra',' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
disp(strcat('se = ', num2str(std(bootstat(:,PostStim),0,'all'))));

%% REVIEWER ANALYSES %%

%% Cell type PSP histograms/variance:

uniqueTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(all_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);

%all_psp_means_std = zeros(size(uniqueTypes,1), 2);
figure('Position', [10 10 1000 900]);
%figure;
for i = 1:size(uniqueTypes)
    if isnan(uniqueTypes{i})
        continue
    end
    ax = subplot(ceil(length(uniqueTypes)/5), 5, i);
    matchType = strcmpi(all_table.classif_data, uniqueTypes{i});
    matchSlice = ismember(all_table.experiment', contra)';
    all_psps = all_response_amps(logical(matchType & matchSlice));
    num_cells = sum(logical(matchType & matchSlice));
    
    cols = hex2rgb(uniqueCols{i});
    
    if num_cells > 1
        histogram(-all_psps, 10, 'FaceColor', cols, 'EdgeColor', 'w');
        title(strcat(uniqueTypes{i}, ', n=', num2str(num_cells)));
        %scatter(mean(all_psps), (std(all_psps)^2), 100, cols, 'filled');
        %hold on;
        all_psp_means_std(i, 1) = mean(all_psps);
        all_psp_means_std(i, 2) = (std(all_psps)^2);
    end
end

%% Replicate PSP variance:

uniqueTypes = unique(all_table.subclass_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(all_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);

figure;
for i = 1:2
    ax = subplot(1, 2, i);
    exp_oi = thalamus{i};

    for j = 1:size(uniqueTypes)
        if isnan(uniqueTypes{j})
            continue
        end
        matchExp = strcmpi(all_table.experiment, exp_oi);
        matchType = strcmpi(all_table.subclass_data, uniqueTypes{j});
        all_psps = all_response_amps(logical(matchType & matchExp));
        num_cells = sum(logical(matchType & matchExp));
        
        cols = hex2rgb(uniqueCols{j});
        
        if num_cells > 1
            scatter(mean(all_psps), (std(all_psps)^2), 100, cols, 'filled');
            hold on;
        end
    end 
end

%% Responder vs non-responder average PSP traces:

figure;

ax = subplot(1, 2, 1);
matchExp = ismember(responders_table.experiment, contra);
cells = responders_table.snip_data(logical(matchExp));
mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
mtrace = cell2mat(mtrace)';
mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
smtrace = -mtrace;
smtrace = smoothdata(smtrace, 'gaussian', 8);

n = sum(~isnan(mean(smtrace,1)),2);
[ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.3); xlim([0 180]);
p.EdgeColor = cols;
hold on;
plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2)))); ylim([-0.5, 6]); hold on;

ax = subplot(1, 2, 2);
matchExp = ismember(non_responders_table.experiment, contra);
cells = non_responders_table.snip_data(logical(matchExp));
mtrace = cellfun(@(x) mean(x, 3), cells, 'UniformOutput', false);
mtrace = cell2mat(mtrace)';
mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
smtrace = -mtrace;
smtrace = smoothdata(smtrace, 'gaussian', 8);

n = sum(~isnan(mean(smtrace,1)),2);
[ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.3); xlim([0 180]);
p.EdgeColor = cols;
hold on;
plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2)))); ylim([-0.5, 6]);


%% ANOVA for exp vs cell type

% uniqueTypes = unique(all_table.classif_data, 'stable');
% type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
% [type_order, idx] = sortrows(type_order, [1 2]);
% 
% boot_ci_means = NaN(numel(uniqueTypes),2);
% 
% for i = 1:size(uniqueTypes)
% 
%     matchType = strcmpi(all_table.classif_data, uniqueTypes{i});
%     matchExp = strcmpi(all_table.experiment, '0802');
%     cells = all_mm_norm_amps(logical(matchType & matchExp));
%     mPSP = mean(cells, 'all', 'omitNan');
%     boot_ci_means(i,2) = mPSP;
% 
% end
% 
% contra_ci_means = boot_ci_means;
% thal_ci_means = boot_ci_means;

exp_one = '0730';
exp_two = '0802';

stat_struct = vertcat(all_mm_norm_amps(strcmpi(all_table.experiment, exp_one)), all_mm_norm_amps(strcmpi(all_table.experiment, exp_two)));
g1 = vertcat(all_table.classif_data(strcmpi(all_table.experiment, exp_one)), all_table.classif_data(strcmpi(all_table.experiment, exp_two)));
g2 = vertcat(all_table.experiment(strcmpi(all_table.experiment, exp_one)), all_table.experiment(strcmpi(all_table.experiment, exp_two)));
[~,~,stats] = anovan(stat_struct, {g1,g2});
[results,~,~,gnames] = multcompare(stats, 'Dimension', [1 2]);

tbl_struct = vertcat(all_mm_norm_amps(strcmpi(all_table.experiment, exp_one)), all_mm_norm_amps(strcmpi(all_table.experiment, exp_two)));
g1 = vertcat(all_table.classif_data(strcmpi(all_table.experiment, exp_one)), all_table.classif_data(strcmpi(all_table.experiment, exp_two)));
g2 = vertcat(all_table.experiment(strcmpi(all_table.experiment, exp_one)), all_table.experiment(strcmpi(all_table.experiment, exp_two)));
tbl = table(tbl_struct, g1, g2, 'VariableNames', {'psp', 'type', 'exp'});
figure; boxchart(categorical(tbl.type), tbl.psp, 'GroupByColor', tbl.exp);

%% Cell type traces, with SEMs, with each experiment normalized to the max mean response:

uniqueTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(all_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);

boot_ci_means = NaN(numel(uniqueTypes),1);
boot_ci_sem = NaN(numel(uniqueTypes),1);
boot_ci_n = NaN(numel(uniqueTypes), 1);
whole_boot_ci = NaN(size(uniqueTypes,1), 50000);

scale_thalamus = [0.5824, 1.3736];
scale_contra = [1.031, 3.1802];

figure('Position', [10 10 1000 900]);

for i = 1:size(uniqueTypes)
    scale_factor = scale_contra;
    if isnan(uniqueTypes{i})
        continue
    end
    ax = subplot(ceil(length(uniqueTypes)/5), 5, i);
    matchType = strcmpi(all_table.classif_data, uniqueTypes{i});
    matchSlice = ismember(all_table.experiment', contra)';
    match_exp_1 = strcmpi(all_table.experiment, contra{1});
    match_exp_2 = strcmpi(all_table.experiment, contra{2});

    cells_1 = all_table.snip_data(logical(matchType & matchSlice & match_exp_1));
    num_cells_1 = sum(logical(matchType & matchSlice & match_exp_1));
    mtrace_1 = cellfun(@(x) mean(x, 3)/scale_factor(1), cells_1, 'UniformOutput', false);
    mtrace_1 = cell2mat(mtrace_1)';
    cells_2 = all_table.snip_data(logical(matchType & matchSlice & match_exp_2));
    num_cells_2 = sum(logical(matchType & matchSlice & match_exp_2));
    mtrace_2 = cellfun(@(x) mean(x, 3)/scale_factor(2), cells_2, 'UniformOutput', false);
    mtrace_2 = cell2mat(mtrace_2)';
    
    mtrace = [mtrace_1 mtrace_2];
    mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
    smtrace = -mtrace;
    smtrace = smoothdata(smtrace, 'gaussian', 8);
    
    if size(smtrace,2) == 0
        continue
    end
    cols = hex2rgb(uniqueCols{i});
    
    if size(smtrace,2) > 1
        n = sum(~isnan(mean(smtrace,1)),2);
        [ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
        upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
        lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
        whole_boot_ci(i,:) = mean(bootstat(:,PostStim), 2, 'omitnan');

        %p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.3); xlim([0 180]);
        %p.EdgeColor = cols;
        %hold on;
        %plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2);  ylim([-0.5, 2]);
        histogram(mean(smtrace(PostStim,:), 1, 'omitnan'), 'BinWidth', 0.4, 'FaceColor', cols, 'EdgeColor', 'w'); xlim([-0.4 4.4])
        title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
        disp(strcat(uniqueTypes{i},' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
        %max_idx = find(mean(smtrace, 2, 'omitnan') == max(mean(smtrace, 2, 'omitnan')));
        disp(strcat('se = ', num2str(std(bootstat(:,PostStim),0,'all'))));
        boot_ci_means(i) = mean(smtrace(PostStim,:), 'all', 'omitnan');
        boot_ci_sem(i) = std(bootstat(:,PostStim),0,'all');
        boot_ci_n(i) = size(smtrace, 2);
        yline(0, '--');
    else
        plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-0.5, 2]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
        disp(strcat(uniqueTypes{i},' = ',num2str(max(smtrace))));
        yline(0, '--');
    end
end

%% %% Response traces, with SEMs, with each experiment normalized to the max mean response:

edges = [-1, -0.9, -0.7, -0.4, -0.25, -0.075, 0]; % Laminar bins based on subjective boundaries. First bin is L6b, last bin is L1
binned_depth_all = discretize(all_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'L5b', 'L4/5a', 'L2/3', 'L1'});
binned_depth_pos = discretize(responders_table.depth_data, edges, 'categorical', {'L6b', 'L6a', 'L5b', 'L4/5a', 'L2/3', 'L1'});
all_depths = flip(unique(binned_depth_all));

scale_thalamus = [0.3749, 0.45869];
scale_contra = [0.66775, 2.2008];

figure('Position', [10 10 1000 200]);
count = 1;
for i = 1:(numel(all_depths))
    scale_factor = scale_thalamus;
    ax1 = subplot(1, numel(all_depths), count);
    matchSlice = ismember(all_table.experiment', thalamus)';
    matchDepth = binned_depth_all == all_depths(i);
    matchType = strcmpi(all_table.subclass_data, '052 Pvalb Gaba'); %053 Sst Gaba

    match_exp_1 = strcmpi(all_table.experiment, thalamus{1});
    match_exp_2 = strcmpi(all_table.experiment, thalamus{2});

    cells_1 = all_table.snip_data(logical(matchDepth & matchSlice & matchType & match_exp_1));
    cells_2 = all_table.snip_data(logical(matchDepth & matchSlice & matchType & match_exp_2));
    num_cells_1 = sum(logical(matchDepth & matchSlice & match_exp_1));
    num_cells_2 = sum(logical(matchDepth & matchSlice & match_exp_2));
    if ~(num_cells_1 + num_cells_2)
        count = count + 1;
        continue
    else
        mtrace_1 = cellfun(@(x) mean(x, 3)/scale_factor(1), cells_1, 'UniformOutput', false);
        mtrace_1 = cell2mat(mtrace_1)';
        mtrace_2 = cellfun(@(x) mean(x, 3)/scale_factor(2), cells_2, 'UniformOutput', false);
        mtrace_2 = cell2mat(mtrace_2)';

        mtrace = [mtrace_1 mtrace_2];
        mtrace(Stim,:) = zeros(numel(Stim),size(mtrace,2));
        smtrace = -mtrace;
        smtrace = smoothdata(smtrace, 'gaussian', 8);
        if size(smtrace,2) == 0
            continue
        end
    
        cols = 'b';
        
        if size(smtrace,2) > 1
            n = sum(~isnan(mean(smtrace,1)),2);
            [ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
            upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
            lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
            %upperci = smoothdata(ci(2,:)', 'gaussian', 8);
            %lowerci = smoothdata(ci(1,:)', 'gaussian', 8);
            %upperci = mean(smtrace, 2, 'omitnan')+tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
            %lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.90,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
            plot(upperci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]);
            hold on;
            plot(lowerci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-0.5, 3]);
            hold on;
            p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.5);
            p.EdgeColor = cols;
            hold on;
            plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(string(all_depths(i)), ', n=', num2str(size(smtrace, 2))));
            disp(strcat(uniqueTypes{i},' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
            %max_idx = find(mean(smtrace, 2, 'omitnan') == max(mean(smtrace, 2, 'omitnan')));
            disp(strcat('se = ', num2str(std(bootstat(:,PostStim),0,'all'))));
            %plot(smtrace, 'color', cols, 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 5]); title(strcat(excit_types{i}, ', n=', num2str(size(smtrace, 2)), {sprintf('\n')}, string(all_depths(k))));
            %hold on;
            %plot(mean(smtrace, 2, 'omitnan'), 'color', 'black', 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 4]);

        else
            plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-0.5, 3]); title(strcat(string(all_depths(i)), ', n=', num2str(size(smtrace, 2))));
        end
        count = count + 1;
    end
end

%% Plotting optimal trial number by cell type (and circuit)

uniqueTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueTypes = uniqueTypes(idx);
responderTypes = categorical(responders_table.classif_data, uniqueTypes);
circuits = categorical(string(ismember(responders_table.experiment, contra)));

cells = logical(ismember(responders_table.experiment, contra));
trial_nums = cellfun(@(x) size(x, 3), responders_table.snip_data, 'UniformOutput', false);
trial_nums = cell2mat(trial_nums);

plot_table = table(trial_nums(cells), responderTypes(cells), circuits(cells));

% figure; boxchart(plot_table.Var2, plot_table.Var1, 'GroupByColor', plot_table.Var3);
% figure; violinplot(plot_table.Var2, plot_table.Var1, GroupByColor=circuits(cells));
figure; swarmchart(plot_table.Var2, plot_table.Var1, 100, '.');


% Two way anova to explain variance using a linear model:

T = table(trial_nums, responderTypes, circuits, ...
    'VariableNames', {'trial_num','cellType','circuit'});

lm = fitlm(T, 'trial_num ~ cellType + circuit');
anovaTbl = anova(lm, 'components');
disp(anovaTbl)

SS = anovaTbl.SumSq;
terms = anovaTbl.Properties.RowNames;

SS_total = sum(SS);

eta2 = SS / SS_total;

etaTbl = table(terms, eta2, ...
    'VariableNames', {'Term', 'EtaSquared'});
disp(etaTbl)

SS_error = SS(strcmp(terms, 'Error'));

partial_eta2 = SS ./ (SS + SS_error);

partialEtaTbl = table(terms, partial_eta2, ...
    'VariableNames', {'Term', 'PartialEtaSquared'});
disp(partialEtaTbl)

c = unique(circuits);

CircuitMeans = zeros(numel(c), 1);
tn1 = trial_nums;
for i = 1:numel(c)
    CircuitMeans(i) = mean(trial_nums(circuits==c(i)));
    tn1(circuits==c(i)) = tn1(circuits==c(i))-CircuitMeans(i);
end


c = unique(responderTypes);

CircuitMeans = zeros(numel(c), 1);
tn2 = trial_nums;
for i = 1:numel(c)
    CircuitMeans(i) = mean(trial_nums(responderTypes==c(i)));
    tn2(responderTypes==c(i)) = tn2(responderTypes==c(i))-CircuitMeans(i);
end


% Variance in response amplitude explained by cell type and depth

T = table(all_mm_norm_amps(ismember(all_table.experiment, thalamus)), ...
    categorical(all_table.classif_data(ismember(all_table.experiment, thalamus))), ...
    binned_depth_all(ismember(all_table.experiment, thalamus)), ...
    'VariableNames', {'response','cellType', 'depth'});

lm = fitlm(T, 'response ~ cellType + depth');
anovaTbl = anova(lm, 'components');
disp(anovaTbl)

SS = anovaTbl.SumSq;
terms = anovaTbl.Properties.RowNames;

SS_total = sum(SS);

eta2 = SS / SS_total;

etaTbl = table(terms, eta2, ...
    'VariableNames', {'Term', 'EtaSquared'});
disp(etaTbl)