function distances = load_stim_distances(param,fun)
    if nargin<2; fun = @mean; end
    fullfnDMD = getFileName(param, 'DMDPatternInfo');
    fullfnRoi = getFileName(param, 'ROI');
    if exist(fullfnDMD,'file')~=2; return; end
    dmdInfo = load(fullfnDMD);

    region = FrameReader(fullfile(param.pth,param.imtoken)).region;
    if exist(getFileName(param, 'OptoData'),'file')~=2; [~,~] = load_stim_series(param); end
    camROI = load(fullfnRoi).roi;
    dmdRoi = extractROIsFromPatterns(param, dmdInfo);

    nDMDRoi = size(dmdRoi, 3);
    nCamROI = size(camROI,3);
    fullCamROIs = false(3200, 3200, nCamROI);
    camDMDPattern = false(3200, 3200, nDMDRoi);
    
    % Generate camDMDPattern
    for j = 1:nDMDRoi
        camDMDPattern(:,:,j) = imwarp(dmdRoi(:,:,j), invert(dmdInfo.excitationToStimRegistration.T2), 'OutputView', imref2d([3200, 3200])) > 0.5;
        camDMDPattern(:,:,j) = camDMDPattern(:,:,j)';
    end

    for j = 1:nCamROI
        fullCamROIs(region(3):region(4),region(1):region(2),j) = camROI(:,:,j);
    end
    
    % Preallocate distance matrix
    distances = zeros(size(camROI, 3), nDMDRoi);
    % Compute the shortest distance between camROI and camDMDPattern
    for i = 1:nCamROI
        % Compute the distance transform for the current camROI
        distanceTransform = bwdist(fullCamROIs(:,:,i));
        
        for j = 1:nDMDRoi
            % Get the current pattern
            patternJ = camDMDPattern(:,:,j);
            
            % Find the minimum distance using the distance transform
            distances(i, j) = fun(distanceTransform(patternJ));
        end
    end
    % j = 4; i=2; figure;imagesc(camDMDPattern(:,:,j) + fullCamROIs(:,:,i)); xlim([1500,2000]); ylim([600,1300]); distances(i,j);
end