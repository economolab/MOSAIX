cd '/Volumes/SynMap Data/2024-08-11 SynMap Contra MCtx/Registered';
responders = readtable('HCR_classification.csv');
responders = table2struct(responders, 'ToScalar', true);

file_list = dir('*Slice*');
file_names = string({file_list.name}');

opto_data = zeros(0,1);
roi_data = zeros(0,1);
classif_data = zeros(0,1);
subclass_data = zeros(0,1);
classif_cols = zeros(0,1);
slice_depth = zeros(0,1);
slice_data = zeros(0,1);
depth_data = zeros(0,1);
cheriff_data = zeros(0,1);
L1_cheriff = zeros(0,1);
soma_cheriff = zeros(0,1);
cheriff_abs = zeros(0,1);
cheriff_pos = zeros(0,1);
for i = 1:size(file_names, 1)

    sub_responders = find(responders.slice == i);
    sub_ids = responders.cell_ids(sub_responders);
    sub_classifs = string(responders.supertype);
    sub_subclass = string(responders.subclass);
    sub_classifs = sub_classifs(sub_responders);
    sub_subclass = sub_subclass(sub_responders);
    sub_cols = string(responders.classif_cols);
    sub_cols = sub_cols(sub_responders);
    sub_slice = responders.slice(sub_responders);
    sub_slice_depth = responders.slice_depth(sub_responders);
    sub_depth = responders.depth(sub_responders);

    fn = fullfile(pwd, file_names(i,:), 'TS files');
    ts_files = dir(fn);
    ts_files = string({ts_files.name}');
    
    tss = regexp(ts_files, '^TS');
    ts_idx = find(~cellfun(@isempty,tss));
    TSData1 = load(fullfile(fn, ts_files(ts_idx(1))));
    TSData2 = load(fullfile(fn, ts_files(ts_idx(2))));
    if i == 1
        ts_data = cell2mat(vertcat(TSData1.ts, TSData2.ts));
        %ts_data = cell2mat(TSData.ts);
    else
        new_ts = cell2mat(vertcat(TSData1.ts, TSData2.ts));
        ts_data = horzcat(ts_data, new_ts);
    end
    size(ts_data, 2)
    
    optos = regexp(ts_files, '^Opto');
    opto_idx = find(~cellfun(@isempty,optos));
    stimData1 = load(fullfile(fn, ts_files(opto_idx(1))));
    stimData2 = load(fullfile(fn, ts_files(opto_idx(2))));
    opto_data = vertcat(opto_data, find(stimData1.stim == 1));
    opto_data = vertcat(opto_data, find(stimData2.stim == 1) + size(stimData1.stim,1));

    direct_cher = regexp(ts_files, 'CheRiff_composite_table');
    direct_idx = find(~cellfun(@isempty,direct_cher));
    direct_anno = readtable(fullfile(fn, ts_files(direct_idx(2))));

    rois = regexp(ts_files, '^ROI');
    roi_idx = find(~cellfun(@isempty,rois));
    ROIData = load(fullfile(fn, ts_files(roi_idx(2))));
    keepIDs = find(ismember(sub_ids, ROIData.roiIndices));
    keepIDs = sub_ids(keepIDs);
   
    roi_data = vertcat(roi_data, ROIData.roiIndices);
    
    display(size(roi_data, 1))

    cheriff_pos = vertcat(cheriff_pos, direct_anno.UserLabel(keepIDs));
    cher = regexp(ts_files, 'CheRiff_Input');
    cher_idx = find(~cellfun(@isempty,cher));
    cheriff = readtable(fullfile(fn, ts_files(cher_idx(2))));
    cheriff_data = vertcat(cheriff_data, cheriff.normAllLayers(keepIDs));
    cheriff_abs = vertcat(cheriff_abs, cheriff.allLayers(keepIDs));
    L1_cheriff = vertcat(L1_cheriff, cheriff.normL1(keepIDs));
    soma_cheriff = vertcat(soma_cheriff, cheriff.normAtSoma(keepIDs));

    classif_data = vertcat(classif_data, sub_classifs(keepIDs));
    subclass_data = vertcat(subclass_data, sub_subclass(keepIDs)); 
    classif_cols = vertcat(classif_cols, sub_cols(keepIDs));
    slice_data = vertcat(slice_data, sub_slice(keepIDs));
    slice_depth = vertcat(slice_depth, sub_slice_depth(keepIDs));
    depth_data = vertcat(depth_data, sub_depth(keepIDs));

end


screenOptoix = zeros(size(ts_data, 1)/400,1);
optoix = unique(opto_data);
stim = 0;
j = 1;
for i = 1:numel(optoix)
    if optoix(i) > stim + 300
        screenOptoix(j) = optoix(i);
        stim = optoix(i);
        j = j + 1;
    else
        continue
    end
end
optoix = nonzeros(screenOptoix);

Nroi = numel(roi_data);

% respondersTrim = responders;
% respondersTrim.cell_ids = responders.cell_ids(ismember(responders.cell_ids, rois));
% respondersTrim.classif = responders.classif(ismember(responders.cell_ids, rois));
% respondersTrim.classif_cols = responders.classif_cols(ismember(responders.cell_ids, rois));

% cheriffToKeep = table2struct(cheriff(:,1:3), 'ToScalar',true);
% cheriffToKeep.cellID = cheriff.cellID(ismember(cheriff.cellID, rois));
% cheriffToKeep.allLayers = cheriff.allLayers(ismember(cheriff.cellID, rois));
% cheriffToKeep.normAllLayers = cheriff.normAllLayers(ismember(cheriff.cellID, rois));

tsAll = ts_data;
F = zeros(size(tsAll));
for i = 1:Nroi
    F(:,i) = -(tsAll(:,i) - median(tsAll(:,i)))./-(median(tsAll(:,i))); % df/f
    F(:,i) = F(:,i) - movmedian(F(:,i), 1000); % subtract bleaching and instability
end

df = 100.*F;

tinterp = -35:1:200;
PreStim = 1:30;
%PostStim = 25:65;
%PostEPSP = 150:180;
%LongPostEPSP = 130:180;


disnip = zeros(numel(tinterp), numel(optoix), Nroi); % For opto stim ix

for i = 1:numel(optoix)
   
    tstim = optoix(i);
    for j = 1:Nroi
        %df = tsAll(:,j);
        dimed = mean(df(tstim+tinterp([PreStim]), j));
        %disnip(:, i, j) = interp1(1:numel(df(:,j)), df(:,j), tstim+tinterp, 'linear');
        disnip(:, i, j) = df(tstim+tinterp,j)-dimed;%-dimed; %ts{j}(:,j)
    end
end

disnip(:,99:100,:) = [];

%response_amps = squeeze(-mean(disnip(PostStim,:,:),[1 2]));
% figure; 
% scatter(slice_data, response_amps); xlim([2 8]); hold on;
% scatter(slice_data, cheriff_abs/200);
% figure; scatter(cheriff_abs, response_amps, '*'); ylim([-0.5 1.5]);

uniqueTypes = unique(classif_data, 'stable');
[uniqueTypes, idx] = sort(uniqueTypes);
uniqueCols = unique(classif_cols, 'stable');
uniqueCols = uniqueCols(idx);

%% Responses plotted and normalized by CheRiff signal in the column %%
% figure;
% 
% for i = 1:length(uniqueTypes)
%     if isnan(uniqueTypes{i})
%         continue
%     end
%     ax = subplot(ceil(length(uniqueTypes)/5), 5, i);
%     matchType = strcmpi(classif_data, uniqueTypes{i});
%     direct_cher = ~cheriff_pos;
%     matchSlice = (ismember(slice_data, [2,3]));
%     cells = disnip(:,:,logical(direct_cher & matchType & matchSlice));
%     typeData = cells(:,1:(ceil(size(disnip,2)/1)),:); % & responders.ResponseSNR > 0
%     cheriffForCells = cheriff_data(logical(direct_cher & matchType & matchSlice));
%     mtrace = squeeze(mean(typeData, 2));
%     smtrace = -mtrace;
%     %smtrace = smtrace./cheriffForCells';
%     smtrace = smoothdata(smtrace, 'gaussian', 4);
%     snr = mean(mtrace(PostStim,:), 1) < -(std(mtrace(PostEPSP,:),0,1)/2); %*1.645
%      % 
%      % for j = 1:size(smtrace, 2)
%      %     if max(smtrace(PostStim,j)) > 5 %|| ~snr(j)
%      %         smtrace(:,j) = NaN(size(smtrace, 1), 1);
%      %     end
%      % end
%     n = sum(~isnan(mean(smtrace,1)),2);
%     upperci = mean(smtrace, 2, 'omitnan')+tinv(0.95,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
%     lowerci = mean(smtrace, 2, 'omitnan')-tinv(0.95,n-1)*std(smtrace,0,2, 'omitmissing')/sqrt(n);
% 
%     % cols = 'black';
%     % if ~strcmp(uniqueTypes{i}, 'ZNo markers')
%     %     cols = hex2rgb(uniqueCols{i});
%     % end
%     cols = hex2rgb(uniqueCols{i});
% 
%     %plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-1.0, 2.5]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(typeData, 3)))); hold on;
%     if size(smtrace,2) > 1
%         %plot(upperci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-1, 1]);
%         %hold on;
%         %plot(lowerci, 'color', cols, 'LineWidth', 0.1); xlim([0 180]); ylim([-1, 1]);
%         %hold on;
%         %p = fill([(1:size(lowerci,1)).'; flip(1:size(upperci,1)).'],[lowerci; flip(upperci)], cols, 'FaceAlpha', 0.5);
%         %p.EdgeColor = cols;
%         %hold on;
%         %plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); title(strcat(uniqueTypes{i}, ', n=', num2str(size(typeData, 3))));
%          % 
%         %imagesc(smtrace'); clim([-1.0 5]); colormap(my_cmap); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
% 
%         plot(smtrace, 'color', cols, 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 3]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
%         hold on;
%         plot(mean(smtrace, 2, 'omitnan'), 'color', 'black', 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 3]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(typeData, 3))));
% 
%         %plot(mean(smtrace, 2, 'omitnan'), 'color', cols, 'LineWidth', 2); xlim([0 180]); ylim([-1.0, 2.0]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(typeData, 3))));
% 
%     else
%         plot(smtrace, 'color', cols, 'LineWidth', 1.5); xlim([0 180]); ylim([-1, 3]); title(strcat(uniqueTypes{i}, ', n=', num2str(size(smtrace, 2))));
%     end
% end
% 
% %% Responses plotted by depth: %%
% 
% matchSlice = (slice_data == 3 | slice_data == 4 | slice_data == 6 | slice_data == 7);
% matchCell = ~cheriff_pos; %ismember(classif_data, classif_data);  %(classif_data == '0115 L6 CT CTX Glut_2');
% [sorted_depth, I] = sort(unique(round(depth_data(logical(matchSlice & matchCell)), 2)), 'descend');
% [lines, ~, subs] = unique(round(depth_data(logical(matchSlice & matchCell)), 2));
% sorted_unique_depth = lines(I);
% response_depth_mean = accumarray(subs, response_amps(logical(matchSlice & matchCell)), [], @sum);
% response_depth_mean = response_depth_mean(I);
% 
% figure; % This will show response amps by depth, normalized for mean column cheriff
% res_plot = scatter(response_amps(logical(matchSlice & matchCell))./cheriff_data(logical(matchSlice & matchCell)), ...
%     depth_data(logical(matchSlice & matchCell)), 65, ...
%     hex2rgb(classif_cols(logical(matchSlice & matchCell))), 'filled'); 
% hold on;
% ylim([-1.0 0]); 
% xlim([-1.0 6.0]);
% xline(median(response_amps(logical(matchSlice & matchCell))));
% res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Slice', slice_data(logical(matchSlice & matchCell)));
% res_plot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Roi Num', roi_data(logical(matchSlice & matchCell)));
% 
% figure; imagesc(ones(size(sorted_unique_depth)), -sorted_unique_depth, response_depth_mean); ylim([0, 1.0]); clim([0 10]); colormap 'hot';
% 
% figure; scatter(slice_data(matchCell), response_amps(matchCell), 65, hex2rgb(classif_cols(matchCell)), 'filled'); xlim([2 8]);

%% Save the disnip object with all associated metadata

experiment = cell(size(roi_data));
experiment(:) = {'0811'};
snip_data = permute(disnip, [3 1 2]);
snip_data = num2cell(snip_data,[2 3]);
to_write = table(snip_data, roi_data, classif_data, subclass_data, classif_cols, ...
    cheriff_abs, cheriff_data, L1_cheriff, soma_cheriff, cheriff_pos, ...
    depth_data, slice_data, slice_depth, experiment);
save('Responses_by_type_0811.mat', 'to_write')
