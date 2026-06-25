function rois = extractROIsFromPatterns(param,dmdInfo)
if exist("dmdInfo","var")
    dmdPatterns = dmdInfo.dmdPatterns;
else
    load(getFileName(param,"DMDPatternInfo"),"dmdPatterns");
end

% Preallocate the output array
labeled_slices = zeros(size(dmdPatterns));

currTotal  = 0;
% Loop through each slice
for k = 1:size(dmdPatterns, 3)
    % Label the connected components in the k-th slice
    [slice, num] = bwlabel(dmdPatterns(:, :, k), 4); % Using 4-connectivity
    % num returns the number of unique connected components found
    slice(slice~=0) = slice(slice~=0) + currTotal;
    labeled_slices(:, :, k) = slice;
    
    currTotal = currTotal+num;
end

imageToGet = sum(labeled_slices,3)>0;
currImage = zeros(size(imageToGet));

[counts,~] = histcounts(labeled_slices,currTotal+1);
counts = counts(2:end);
roiIndices = [];
i = 1;
imageToGet(:,500:end) = 0;
while any((currImage>0)~=(imageToGet>0),'all')
    [x,y] = find(imageToGet & (~currImage),1);
    roiVals = squeeze(labeled_slices(x,y,:));
    roiVals(roiVals==0)=[];
    [~,dex] = min(counts(roiVals));
    roiVal = roiVals(dex);
    roiIndices(end+1) = roiVal;
    currImage(sum(labeled_slices==roiVal,3)>0) = roiVal;
    rois(:,:,i) = sum(labeled_slices==roiVal,3)>0;
    i = i+1;
end