%% Load data and define parameters

% Run from the repo root (like the other MATLAB scripts); paths below are relative to it.
paths = jsondecode(fileread('config/paths.json'));
output_dir = fullfile(pwd, paths.output_dir);

% Load analysis output tables
load(fullfile(output_dir, 'all_table.mat'), 'all_table', 'all_med_norm_amps');
load(fullfile(output_dir, 'responders_table.mat'), 'responders_table');
load(fullfile(output_dir, 'non_responders_table.mat'), 'non_responders_table');

% Reconstruct full response_table (all cells before responder filtering)
response_table = all_table;
response_table = removevars(response_table, 'is_responder');

% Experiment groupings
thalamus = {'0730', '0802'};
contra = {'0811', '0814'};

% Stimulus timing windows
PreStim = 1:32;
Stim = 34:39;
PostStim = 43:63;
PostEPSP = 175:200;

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
    %Responders
    ax = subplot(2, numel(unique_exps), i);
    exp_idx = strcmpi(responders_table.experiment, unique_exps{i});
    exp_cells = responders_table.snip_data(logical(exp_idx));
    mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
    max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
    max_rows = cell2mat(max_rows);
    norm_rows = cellfun(@(x) x./abs(min(x(PostStim))), mean_rows, 'UniformOutput', false);
    norm_rows = cell2mat(norm_rows);
    %[sort_arr, sort_idx] = sort(max_rows);
    sorted_types = sort(responderTypes(logical(exp_idx)));
    [sort_arr, sort_idx] = sortrows([[grp2idx(responderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
    imagesc(-norm_rows(sort_idx,:)); clim ([-1 1]); colormap(cmi); title(strcat('Responders, n = ', num2str(sum(exp_idx))));
    
    %Non-responders
    ax = subplot(2, numel(unique_exps), i+numel(unique_exps));
    exp_idx = strcmpi(non_responders_table.experiment, unique_exps{i});
    exp_cells = non_responders_table.snip_data(logical(exp_idx));
    mean_rows = cellfun(@(x) mean(x, 3), exp_cells, 'UniformOutput', false);
    max_rows = cellfun(@(x) min(x(PostStim)), mean_rows, 'UniformOutput', false);
    max_rows = cell2mat(max_rows);
    norm_rows = cellfun(@(x) x./abs(max(x(PostStim))), mean_rows, 'UniformOutput', false);
    norm_rows = cell2mat(norm_rows);
    %[sort_arr, sort_idx] = sort(max_rows);
    sorted_types = sort(nonresponderTypes(logical(exp_idx)));
    [sort_arr, sort_idx] = sortrows([[grp2idx(nonresponderTypes(logical(exp_idx)))] [max_rows]], [1 2]);
    imagesc(-norm_rows(sort_idx,:)); clim ([-1 1]); colormap(cmi); title(strcat('Non-responders, n = ', num2str(sum(exp_idx))));

end
colorbar;

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

cols = 'r';

n = sum(~isnan(mean(smtrace,1)),2);
[ci, bootstat] = bootci(50000, {@(x) mean(x,1,'omitnan'), -mtrace'}, Options=statset(UseParallel=true));
upperci = mean(smtrace, 2, 'omitnan') + smoothdata(std(bootstat,0,1)', 'gaussian', 8);
lowerci = mean(smtrace, 2, 'omitnan') - smoothdata(std(bootstat,0,1)', 'gaussian', 8);
p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.3); xlim([0 180]);
p.EdgeColor = cols;
hold on;
plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); ylim([-0.5, 6]); hold on;

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
plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); ylim([-0.5, 6]);


%% Response amplitudes plotted by laminar depth %%

plotting_table = all_table(logical(ismember(all_table.experiment, contra)'),:);
figure; % This will show response amps by depth, normalized for mean column cheriff
res_plot = scatter(log2(all_med_norm_amps(logical(ismember(all_table.experiment, contra)'))), plotting_table.depth_data, 65, hex2rgb(plotting_table.classif_cols), 'filled'); 
hold on;
ylim([-1.0 0]); 
xlim([-0.5 8]);
xline([1.0]);
title(strcat('n =  ', num2str(size(plotting_table,1)), ' cells'))
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Exp', plotting_table.experiment);
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Slice', plotting_table.slice_data);
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Roi Num', plotting_table.roi_data);
res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('CellType', plotting_table.classif_data);

%% Stacked bar of responders vs non responders: %%

allTypes = unique(all_table.classif_data, 'stable');
allCols = unique(all_table.classif_cols, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
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
b = bar(percents(:,1:3)', 'stacked', 'FaceColor','flat'); xticklabels({'All cells','Non Responders','Responders'});
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

%% Response amplitudes relative to columnar cheriff %%

allTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(allTypes, ' ')] [extractBefore(allTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
allCols = unique(responders_table.classif_cols, 'stable');
allTypes = allTypes(idx);
grouping = categorical(all_table.classif_data, allTypes);
color_grouping = ismember(all_table.experiment, contra);
combined_color_grouping = cat(1,grouping(color_grouping), grouping(color_grouping));
combined_amps = cat(1, all_med_norm_amps(color_grouping), all_med_norm_amps(color_grouping)./all_table.cheriff_data(color_grouping));
cheriff_group = cat(2, repelem(1, numel(combined_amps)/2), repelem(2, numel(combined_amps)/2));

figure('Position', [10 10 1000 900]);
ax = subplot(2, 1, 1);
boxchart(combined_color_grouping, combined_amps, 'GroupByColor', cheriff_group);
title('Responses'); ylabel('Median-normalized response amp');
%ylim([-1 12]);

%% %% Response traces across layers, with SEMs, with each experiment normalized to the max mean response:

uniqueTypes = unique(all_table.classif_data, 'stable');
type_order = [[extractAfter(uniqueTypes, ' ')] [extractBefore(uniqueTypes, ' ')]];
[type_order, idx] = sortrows(type_order, [1 2]);
uniqueCols = unique(all_table.classif_cols, 'stable');
uniqueTypes = uniqueTypes(idx);
uniqueCols = uniqueCols(idx);

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
    matchType = strcmpi(all_table.subclass_data, '052 Pvalb Gaba'); % Can be excluded as a filter

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
            plot(upperci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]);
            hold on;
            plot(lowerci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-0.5, 3]);
            hold on;
            p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.5);
            p.EdgeColor = cols;
            hold on;
            plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(string(all_depths(i)), ', n=', num2str(size(smtrace, 2))));

        else
            plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-0.5, 3]); title(strcat(string(all_depths(i)), ', n=', num2str(size(smtrace, 2))));
        end
        count = count + 1;
    end
end

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

        p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.3); xlim([0 180]);
        p.EdgeColor = cols;
        hold on;
        plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2);  ylim([-0.5, 2]);
        title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
        disp(strcat(uniqueTypes{i},' = ',num2str(mean(smtrace(PostStim,:), 'all', 'omitnan'))));
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
