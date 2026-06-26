function [ts_out, opto_out, roi_out, keepIDs, cher_data, keepIDs_idx] = load_slice_data(fn, ts_files, slice_config, output_dir, slice_label)
% LOAD_SLICE_DATA Load time series, optogenetic, ROI, and CheRiff data for one slice.
%
%   [ts_out, opto_out, roi_out, keepIDs, cher_data, keepIDs_idx] = load_slice_data(fn, ts_files, slice_config, output_dir, slice_label)
%
%   Inputs:
%     fn           - Path to the slice's 'TS files' directory
%     ts_files     - String array of filenames in fn
%     slice_config - Struct with fields:
%                      ts_idx_list   - which TS file indices to load (e.g., [1,2])
%                      opto_idx_list - which Opto file indices to load
%                      roi_idx_list  - which ROI file indices to load
%                      ts_concat     - 'horzcat_separate', 'vertcat_then_cellmat', or 'single'
%                      sub_ids       - cell IDs from HCR classification for this slice
%     output_dir   - Path to output/<experiment_id>/ where CheRiff CSVs live
%     slice_label  - Slice label string (e.g., 's05') for finding CheRiff files
%
%   Outputs:
%     ts_out      - Time series matrix (timepoints x cells)
%     opto_out    - Vector of stimulus onset indices
%     roi_out     - Vector of ROI indices
%     keepIDs     - Cell IDs that matched ROI data
%     cher_data   - Struct with fields: cheriff_data, cheriff_abs, L1_cheriff,
%                   soma_cheriff, cheriff_pos
%     keepIDs_idx - Indices into sub_ids for each kept cell (preserves
%                   duplicates across ROI files)

    % --- Load time series ---
    tss = regexp(ts_files, '^TS');
    ts_idx = find(~cellfun(@isempty, tss));

    if strcmp(slice_config.ts_concat, 'horzcat_separate')
        % Load each TS file separately and horzcat
        ts_out = [];
        for k = 1:length(slice_config.ts_idx_list)
            TSData = load(fullfile(fn, ts_files(ts_idx(slice_config.ts_idx_list(k)))));
            if iscell(TSData.ts)
                ts_out = horzcat(ts_out, cell2mat(TSData.ts));
            else
                ts_out = horzcat(ts_out, TSData.ts);
            end
        end

    elseif strcmp(slice_config.ts_concat, 'vertcat_then_cellmat')
        % Vertcat cell arrays first, then cell2mat
        all_ts = {};
        for k = 1:length(slice_config.ts_idx_list)
            TSData = load(fullfile(fn, ts_files(ts_idx(slice_config.ts_idx_list(k)))));
            if iscell(TSData.ts)
                all_ts = vertcat(all_ts, TSData.ts);
            else
                all_ts = vertcat(all_ts, {TSData.ts});
            end
        end
        ts_out = cell2mat(all_ts);

    elseif strcmp(slice_config.ts_concat, 'single')
        % Load a single TS file
        TSData = load(fullfile(fn, ts_files(ts_idx(slice_config.ts_idx_list(1)))));
        if iscell(TSData.ts)
            ts_out = cell2mat(TSData.ts);
        else
            ts_out = TSData.ts;
        end
    end

    % --- Load optogenetic stimulus data ---
    optos = regexp(ts_files, '^Opto');
    opto_idx = find(~cellfun(@isempty, optos));

    opto_out = [];
    prev_len = 0;
    for k = 1:length(slice_config.opto_idx_list)
        stimData = load(fullfile(fn, ts_files(opto_idx(slice_config.opto_idx_list(k)))));
        opto_out = vertcat(opto_out, find(stimData.stim == 1) + prev_len);
        if k < length(slice_config.opto_idx_list) && ...
           strcmp(slice_config.ts_concat, 'vertcat_then_cellmat')
            prev_len = size(stimData.stim, 1);
        end
    end

    % --- Load ROI data ---
    rois = regexp(ts_files, '^ROI');
    roi_idx_all = find(~cellfun(@isempty, rois));

    roi_out = [];
    keepIDs = [];
    keepIDs_idx = [];
    for k = 1:length(slice_config.roi_idx_list)
        ROIData = load(fullfile(fn, ts_files(roi_idx_all(slice_config.roi_idx_list(k)))));
        roi_indices = ROIData.roiIndices;
        if size(roi_indices, 2) > size(roi_indices, 1)
            roi_indices = roi_indices';
        end
        roi_out = vertcat(roi_out, roi_indices);

        matched = find(ismember(slice_config.sub_ids, ROIData.roiIndices));
        keepIDs = vertcat(keepIDs, slice_config.sub_ids(matched));
        keepIDs_idx = vertcat(keepIDs_idx, matched);
    end

    % --- Load CheRiff data from output directory ---
    % CheRiff composite table (direct stimulation annotation) from TS files dir
    direct_cher = regexp(ts_files, 'CheRiff_composite_table');
    direct_idx = find(~cellfun(@isempty, direct_cher));
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    direct_anno = readtable(fullfile(fn, ts_files(direct_idx(end))));
    cheriff = readtable(fullfile(output_dir, sprintf('%s_CheRiff_Input_to_All_Rois.csv', slice_label)));
    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames');

    cher_data.cheriff_pos = direct_anno.UserLabel(keepIDs);
    cher_data.cheriff_data = cheriff.normAllLayers(keepIDs);
    cher_data.cheriff_abs = cheriff.allLayers(keepIDs);
    cher_data.L1_cheriff = cheriff.normL1(keepIDs);
    cher_data.soma_cheriff = cheriff.normAtSoma(keepIDs);

end
