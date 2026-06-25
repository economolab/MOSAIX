function [correspondence,distances] = load_stim_correspondence(param,useDMDRoi)
    if nargin<2 || isempty(useDMDRoi); useDMDRoi = 0; end
    
    fullfnDMD = getFileName(param, 'DMDPatternInfo');
    if exist(fullfnDMD,'file')~=2; return; end
    dmdInfo = load(fullfnDMD);

    amountToDilate = 0; 

    region = FrameReader(fullfile(param.pth,param.imtoken)).region;
    if exist(getFileName(param, 'OptoData'),'file')~=2; [~,~] = load_stim_series(param); end
    
    if useDMDRoi
        % if isfield(dmdInfo,"roiIndices"); correspondence = dmdInfo.roiIndices{2:end,1}; return; end
        roi = extractROIsFromPatterns(param,dmdInfo);
        fractionRequired = 0.50;
    else
        fullfnRoi = getFileName(param, 'ROI');
        if exist(fullfnRoi,'file')~=2; correspondence=-1; return; end
        camROI = load(fullfnRoi).roi;
        tempCamPatternRoi = zeros(3200,3200);
        outputViewDMD = imref2d(size(zeros(768, 1024)));
        nRoi = size(camROI,3);
        roi = false(768,1024,nRoi);
        for j = 1:nRoi
            tempCamPatternRoi(region(3):region(4),region(1):region(2)) = camROI(:,:,j);
            roi(:,:,j) = imwarp(tempCamPatternRoi',dmdInfo.excitationToStimRegistration.T2, 'OutputView',outputViewDMD)>0;
        end
        fractionRequired = 0.02;
    end
    patterns = dmdInfo.dmdPatterns;
    roi = imdilate(roi, strel('disk', amountToDilate)) > 0;

    
    roiSums = squeeze(sum(sum(roi, [1,2]), 2));
    sums = reshape(patterns,[],size(patterns,3),1)'*reshape(roi,[],size(roi,3),1); %reshape to be [patternID x pix] * [pix x roi] to get [patternID x roi]
    stimCorrespondence =  sums>(roiSums'* fractionRequired);
    stimCorrespondence = stimCorrespondence(2:end,:);
    correspondence = stimCorrespondence;
    if ~useDMDRoi
        save(fullfnRoi,"stimCorrespondence","-append"); 
        optoData = load(getFileName(param, 'OptoData'));
        if ~isfield(optoData,'stimOffset')
            stimOffset = 0;
            save(optoFN,"stimOffset","-append"); 
        end
    end

    if nargout == 2
        distances = load_stim_distances(param);
    else
        distances = [];
    end
end