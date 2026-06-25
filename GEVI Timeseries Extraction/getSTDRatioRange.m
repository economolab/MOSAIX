function [rangeSTDs,stdRatio] = getSTDRatioRange(mov,sigmaInner,sigmaOuter,outerRadius,temporalFilterSigma)
    if ~exist('sigmaInner','var'); sigmaInner = 3;end
    if ~exist('sigmaOuter','var'); sigmaOuter = 3;end
    if ~exist('outerRadius','var'); outerRadius = 25;end
    if ~exist('applyTemporalFilt','var'); temporalFilterSigma = 2;end
    if ~exist('batchSize','var'); batchSize = 500; end

    mov = single(mov);
    movSize = whos('mov');
    g = gpuDevice;

    if movSize.bytes <= g.AvailableMemory/2
        mov = gpuArray(mov);
    end
    fSize = (outerRadius+sigmaOuter)*2+1;
    
    meanMovWithBG = (mean(mov,3));
    stdImAprox = sqrt(meanMovWithBG);
    myFilter = fspecial('gaussian',fSize,sigmaInner)-getGaussianRing(fSize,outerRadius,sigmaOuter);

    frames = size(mov,3);

    for i = 1:batchSize:frames
        idx = i:min(i+batchSize-1, frames);
        mov(:,:,idx) = imfilter(gpuArray(mov(:,:,idx)),gpuArray(myFilter),"circular");
    end

    if temporalFilterSigma>1
        mov = imgaussfilt3(mov,[0.1,0.1,temporalFilterSigma],"FilterSize",[1,1,temporalFilterSigma*2+1]);
    end
    %change to DF/F but use unbackground subtracted for denominator
    mov = (mov-(mean(mov,3)))./meanMovWithBG; 
    
    mov = permute(mov,[3,1,2]); %detrend pixel across time
    mov = detrend(mov,3);
    mov = ipermute(mov,[3,1,2]);

    dffRange = max(mov,[],3)-min(mov,[],3);
    rangeSTDs = dffRange.*stdImAprox; %DFF std is ~ to sqrt(1/intensity)
    rangeSTDs = gather(rangeSTDs);

    stdRatio = std(mov,0,3).*stdImAprox; %DFF std is ~ to sqrt(1/intensity)
    stdRatio = gather(stdRatio);
end