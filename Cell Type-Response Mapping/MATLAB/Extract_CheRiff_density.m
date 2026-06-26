function Extract_CheRiff_density(config_path)
% EXTRACT_CHERIFF_DENSITY  Extract CheRiff laminar signal densities for all slices.
%
%   Extract_CheRiff_density(config_path)
%
%   For each slice in the experiment, loads the CheRiff channel image and
%   Voltron mask, rotates by the config angle, and computes per-ROI signal
%   densities across the full column, L1 (top 80px), and locally around
%   each soma.
%
%   Writes s##_CheRiff_Input_to_All_Rois.csv to output/<experiment_id>/,
%   which is subsequently read by load_slice_data.m.
%
%   Input:
%     config_path - Path to experiment JSON file
%                   (e.g., 'config/experiments/2024-08-14_contra.json')

    % --- Load configuration ---
    config = jsondecode(fileread(config_path));
    paths = jsondecode(fileread('config/paths.json'));
    base_dir = fullfile(paths.mosaix_data_dir, config.data_subdir);
    output_dir = fullfile(pwd, paths.output_dir, config.experiment_id);

    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('Extracting CheRiff densities for: %s\n', config.experiment_label);

    % --- Enumerate slice folders ---
    file_list = dir(fullfile(base_dir, '*Slice*'));
    file_names = string({file_list.name}');
    n_slices = size(file_names, 1);

    fprintf('Found %d slice folders\n', n_slices);

    angles = config.img_lims.angle_deg;
    ymin = config.img_lims.ymin;   % top of cortex (smaller pixel coord)
    ymax = config.img_lims.ymax;   % bottom of cortex (larger pixel coord)

    L1_depth_px = 80;  % L1 boundary: 80 pixels below top of cortex

    for i = 1:n_slices
        slice_dir = fullfile(base_dir, file_names(i));
        fprintf('  Processing %s...\n', file_names(i));

        % --- Find CheRiff channel TIF ---
        cher_files = dir(fullfile(slice_dir, 'MAX_channel_3_*.tif'));
        if isempty(cher_files)
            cher_files = dir(fullfile(slice_dir, 'MAX_channel_*.tif'));
        end
        if isempty(cher_files)
            warning('No CheRiff channel TIF found in %s, skipping.', slice_dir);
            continue;
        end
        tif = tiffreadVolume(fullfile(slice_dir, cher_files(1).name));

        % --- Find Voltron mask TIF ---
        mask_files = dir(fullfile(slice_dir, '*ilastik_masks_AllVoltronPos_consec*'));
        mask_files = mask_files(~startsWith({mask_files.name}, '._'));
        if isempty(mask_files)
            warning('No Voltron mask TIF found in %s, skipping.', slice_dir);
            continue;
        end
        masktif = tiffreadVolume(fullfile(slice_dir, mask_files(1).name));

        % --- Get rotation angle and cortical boundaries ---
        rotAngle = angles(i);
        top_y = ymin(i);              % top of cortex in pixel coords
        bottom_y = ymax(i);           % bottom of cortex in pixel coords
        L1_y = top_y + L1_depth_px;   % lower boundary of L1

        % --- Rotate images ---
        tif = imrotate(tif, -rotAngle, 'nearest', 'crop');
        masktif = imrotate(masktif, -rotAngle, 'nearest', 'crop');

        % --- Extract ROI densities ---
        rois = unique(masktif);
        n_rois = numel(rois) - 1;  % exclude background (0)

        colValues = zeros(n_rois, 1);
        L1Values = zeros(n_rois, 1);
        localValues = zeros(n_rois, 1);

        for j = 2:numel(rois)
            [row, col, ~] = ind2sub(size(masktif), find(masktif == rois(j)));
            centX = floor(mean(col));
            centY = floor(mean(row));

            % 101px-wide strip centered on soma
            leftX = max(centX - 50, 1);
            rightX = min(centX + 50, size(tif, 2));
            upperY = max(centY - 50, 1);
            lowerY = min(centY + 50, size(tif, 1));

            col_pixels = tif(top_y:bottom_y, leftX:rightX);
            L1_pixels = tif(top_y:L1_y, leftX:rightX);
            local_pixels = tif(upperY:lowerY, leftX:rightX);

            colValues(j-1) = mean(col_pixels, 'all');
            L1Values(j-1) = mean(L1_pixels, 'all');
            localValues(j-1) = mean(local_pixels, 'all');
        end

        colNormValues = colValues / max(colValues);
        L1NormValues = L1Values / max(L1Values);
        localNormValues = localValues / max(localValues);

        T = table(double(rois(2:end)), colValues, colNormValues, ...
                  L1Values, L1NormValues, localValues, localNormValues, ...
                  'VariableNames', {'cellID', 'allLayers', 'normAllLayers', ...
                                    'L1', 'normL1', 'atSoma', 'normAtSoma'});

        % --- Determine output filename from folder name ---
        slice_tokens = regexp(file_names(i), 'Slice\s*(\d+)', 'tokens');
        if ~isempty(slice_tokens)
            slice_num_str = sprintf('s%02d', str2double(slice_tokens{1}{1}));
        else
            slice_num_str = sprintf('s%02d', i);
        end

        out_fn = fullfile(output_dir, sprintf('%s_CheRiff_Input_to_All_Rois.csv', slice_num_str));
        writetable(T, out_fn);
        fprintf('    Wrote: %s (%d ROIs)\n', out_fn, n_rois);
    end

    fprintf('Done.\n');
end
