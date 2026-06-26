wd = '/Volumes/SynMap Data/2024-08-14 SynMap Contra MCtx/Registered/Slice 5';
fn = fullfile(wd, 'MAX_channel_3_s05_NT445_JF585_GFP_40X_zstack_bin.tif');
tif = tiffreadVolume(fn);
maskfn = fullfile(wd, 't_s05_ilastik_masks_AllVoltronPos_consec.tiff');
masktif = tiffreadVolume(maskfn);

rotAngle = 12;

figure; imagesc(imrotate(tif, -rotAngle, 'nearest', 'crop'));
maxY = 169;
minY = 1655;
minL1 = 249;
tif = imrotate(tif, -rotAngle, 'nearest', 'crop');
masktif = imrotate(masktif, -rotAngle, 'nearest', 'crop');

rois = unique(masktif);

colValues = zeros(numel(rois)-1,1);
L1Values = zeros(numel(rois)-1,1);
localValues = zeros(numel(rois)-1,1);

for i = 2:numel(rois)
    [row, col, v] = ind2sub(size(masktif), find(masktif == rois(i)));
    %internalRoiValues(i-1) = mean(tif(row, col), 'all');
    centX = floor(mean(col));
    centY = floor(mean(row));
    leftX = centX - 50;
    rightX = centX + 50;
    upperY = centY-50;
    lowerY = centY+50;
    col_pixels = tif(maxY:minY, leftX:rightX);
    L1_pixels = tif(maxY:minL1, leftX:rightX);
    local_pixels = tif(upperY:lowerY, leftX:rightX);
    colValues(i-1) = mean(col_pixels, 'all');
    localValues(i-1) = mean(local_pixels, 'all');
    L1Values(i-1) = mean(L1_pixels, 'all');
end


colNormValues = colValues/max(colValues);
L1NormValues = L1Values/max(L1Values);
localNormValues = localValues/max(localValues);
T = [double(rois(2:numel(rois))), colValues, colNormValues, L1Values, L1NormValues, localValues, localNormValues];
T = array2table(T);
T.Properties.VariableNames(1:7) = {'cellID', 'allLayers', 'normAllLayers', 'L1', 'normL1', 'atSoma', 'normAtSoma'};
writetable(T, fullfile(wd, 's05_CheRiff_Input_to_All_Rois.csv'));
