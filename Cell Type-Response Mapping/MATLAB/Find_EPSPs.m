function Find_EPSPs(config_path)
% FIND_EPSPS  Unified electrophysiology analysis for all MOSAIX experiments.
%
%   Find_EPSPs(config_path)
%
%   Reads an experiment configuration JSON file and extracts synaptic response
%   data paired with HCR cell type classifications.
%
%   Input:
%     config_path - Path to experiment JSON file
%                   (e.g., 'config/experiments/2024-07-30_thalamus.json')
%
%   Output:
%     Saves Responses_by_type_<experiment_id>.mat to the output directory.

    % --- Load configuration ---
    config = jsondecode(fileread(config_path));

    % Load shared paths
    paths = jsondecode(fileread('config/paths.json'));

    base_dir = fullfile(paths.mosaix_data_dir, config.data_subdir);
    output_dir = fullfile(pwd, paths.output_dir, config.experiment_id);
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('Processing experiment: %s\n', config.experiment_label);
    fprintf('Data directory: %s\n', base_dir);

    cd(base_dir);

    % --- Load HCR classification ---
    % Try output directory first, fall back to experiment directory
    classif_path = fullfile(output_dir, 'HCR_classification.csv');
    if ~exist(classif_path, 'file')
        classif_path = 'HCR_classification.csv';
    end
    responders = readtable(classif_path);
    responders = table2struct(responders, 'ToScalar', true);

    % --- Enumerate slices ---
    file_list = dir('*Slice*');
    file_names = string({file_list.name}');
    n_slices = size(file_names, 1);

    fprintf('Found %d slices\n', n_slices);

    % --- Determine per-slice file loading config ---
    % If per_slice_file_counts is defined, use it; otherwise use the default
    if isfield(config.ephys, 'per_slice_file_counts')
        per_slice_configs = config.ephys.per_slice_file_counts;
    else
        % Use the same config for every slice
        default_config = config.ephys.slice_file_counts;
        if iscell(default_config)
            default_config = default_config{1};
        elseif length(default_config) > 1
            default_config = default_config(1);
        end
        per_slice_configs = repmat({default_config}, n_slices, 1);
    end

    % --- Initialize accumulators ---
    ts_data = [];
    opto_data = zeros(0, 1);
    roi_data = zeros(0, 1);
    classif_data = string.empty(0, 1);
    subclass_data = string.empty(0, 1);
    classif_cols = string.empty(0, 1);
    slice_data = zeros(0, 1);
    slice_depth = zeros(0, 1);
    depth_data = zeros(0, 1);
    cheriff_data = zeros(0, 1);
    cheriff_abs = zeros(0, 1);
    soma_cheriff = zeros(0, 1);
    L1_cheriff = zeros(0, 1);
    cheriff_pos = zeros(0, 1);

    % --- Process each slice ---
    for i = 1:n_slices

        slice_num = i + config.slice_offset;
        sub_responders = find(responders.slice == slice_num);
        sub_ids = responders.cell_ids(sub_responders);
        sub_supertype = string(responders.supertype);
        sub_subclass = string(responders.subclass);
        sub_supertype = sub_supertype(sub_responders);
        sub_subclass = sub_subclass(sub_responders);
        sub_cols = string(responders.classif_cols);
        sub_cols = sub_cols(sub_responders);
        sub_slice = responders.slice(sub_responders);
        sub_slice_depth = responders.slice_depth(sub_responders);
        sub_depth = responders.depth(sub_responders);

        fn = fullfile(pwd, file_names(i, :), 'TS files');
        ts_files_list = dir(fn);
        ts_files_str = string({ts_files_list.name}');

        % Get slice-specific config
        if iscell(per_slice_configs)
            sc = per_slice_configs{i};
        else
            sc = per_slice_configs(i);
        end

        % Build the slice_config struct for load_slice_data
        slice_cfg.ts_idx_list = sc.ts_idx(:)';
        slice_cfg.opto_idx_list = sc.opto_idx(:)';
        slice_cfg.roi_idx_list = sc.roi_idx(:)';
        slice_cfg.ts_concat = sc.ts_concat;
        slice_cfg.sub_ids = sub_ids;

        % Extract slice label from folder name (e.g., "Slice 5" -> "s05")
        slice_tokens = regexp(file_names(i), 'Slice\s*(\d+)', 'tokens');
        if ~isempty(slice_tokens)
            slice_label = sprintf('s%02d', str2double(slice_tokens{1}{1}));
        else
            slice_label = sprintf('s%02d', i);
        end

        % Load all data for this slice
        [ts_slice, opto_slice, roi_slice, keepIDs, cher, keepIDs_idx] = ...
            load_slice_data(fn, ts_files_str, slice_cfg, output_dir, slice_label);

        % keepIDs_idx: positions of kept cells within the classification
        % arrays (sub_*), computed per ROI file to preserve duplicates
        % across images.

        opto_data = vertcat(opto_data, opto_slice);
        roi_data = vertcat(roi_data, roi_slice);

        % Accumulate time series
        if isempty(ts_data)
            ts_data = ts_slice;
        else
            ts_data = horzcat(ts_data, ts_slice);
        end

        fprintf('  Slice %d: %d cells, %d total\n', slice_num, size(ts_slice, 2), size(ts_data, 2));

        cheriff_pos = vertcat(cheriff_pos, cher.cheriff_pos);
        cheriff_data = vertcat(cheriff_data, cher.cheriff_data);
        cheriff_abs = vertcat(cheriff_abs, cher.cheriff_abs);
        L1_cheriff = vertcat(L1_cheriff, cher.L1_cheriff);
        soma_cheriff = vertcat(soma_cheriff, cher.soma_cheriff);

        classif_data = vertcat(classif_data, sub_supertype(keepIDs_idx));
        subclass_data = vertcat(subclass_data, sub_subclass(keepIDs_idx));
        classif_cols = vertcat(classif_cols, sub_cols(keepIDs_idx));
        slice_data = vertcat(slice_data, sub_slice(keepIDs_idx));
        slice_depth = vertcat(slice_depth, sub_slice_depth(keepIDs_idx));
        depth_data = vertcat(depth_data, sub_depth(keepIDs_idx));
    end

    % --- Screen opto stimulus indices (keep stimuli >300 frames apart) ---
    screenOptoix = zeros(ceil(size(ts_data, 1) / 400), 1);
    optoix = unique(opto_data);
    stim = 0;
    j = 1;
    for i = 1:numel(optoix)
        if optoix(i) > stim + 300
            screenOptoix(j) = optoix(i);
            stim = optoix(i);
            j = j + 1;
        end
    end
    optoix = nonzeros(screenOptoix);

    Nroi = numel(roi_data);
    fprintf('Total ROIs: %d, Stimulus events: %d\n', Nroi, numel(optoix));

    % --- Compute dF/F with bleaching correction ---
    tsAll = ts_data;
    F = zeros(size(tsAll));
    for i = 1:Nroi
        F(:, i) = -(tsAll(:, i) - median(tsAll(:, i))) ./ -(median(tsAll(:, i)));
        F(:, i) = F(:, i) - movmedian(F(:, i), 1000);
    end
    df = 100 .* F;

    % --- Extract peri-stimulus snippets ---
    tinterp = -35:1:200;
    PreStim = 1:30;

    disnip = zeros(numel(tinterp), numel(optoix), Nroi);
    for i = 1:numel(optoix)
        tstim = optoix(i);
        for j = 1:Nroi
            dimed = mean(df(tstim + tinterp(PreStim), j));
            disnip(:, i, j) = df(tstim + tinterp, j) - dimed;
        end
    end

    % --- Optional: trim disnip columns ---
    if isfield(config.ephys, 'disnip_trim') && ~isempty(config.ephys.disnip_trim)
        trim_cols = config.ephys.disnip_trim(:)';
        disnip(:, trim_cols, :) = [];
        fprintf('Trimmed disnip columns: %s\n', mat2str(trim_cols));
    end

    % --- Save output ---
    experiment = cell(size(roi_data));
    experiment(:) = {config.experiment_id};
    snip_data = permute(disnip, [3 1 2]);
    snip_data = num2cell(snip_data, [2 3]);
    to_write = table(snip_data, roi_data, classif_data, subclass_data, classif_cols, ...
        cheriff_abs, cheriff_data, L1_cheriff, soma_cheriff, cheriff_pos, ...
        depth_data, slice_data, slice_depth, experiment);

    output_file = fullfile(output_dir, sprintf('Responses_by_type_%s.mat', config.experiment_id));
    save(output_file, 'to_write');
    fprintf('Saved: %s\n', output_file);

    % Also save to experiment directory for backwards compatibility
    legacy_file = sprintf('Responses_by_type_%s.mat', config.experiment_id);
    save(legacy_file, 'to_write');
    fprintf('Legacy save: %s\n', fullfile(base_dir, legacy_file));

end
