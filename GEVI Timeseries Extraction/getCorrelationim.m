function correlationim = getCorrelationim(mov,sigmaInner,sigmaOuter,outerRadius,temporalFilter,batchSize)
    if ~exist('sigmaInner','var') || isempty(sigmaInner); sigmaInner = 1;end
    if ~exist('sigmaOuter','var') || isempty(sigmaOuter); sigmaOuter = 3;end
    if ~exist('outerRadius','var') || isempty(outerRadius); outerRadius = 25;end
    if ~exist('applyTemporalFilt','var') || isempty(temporalFilter); temporalFilter = 0;end
    if ~exist('batchSize','var') || isempty(batchSize); batchSize = 500; end

    mov = single(mov);
    fSize = (outerRadius+sigmaOuter)*2+1;
    
    if sigmaOuter>0
        backgroundFilter = getGaussianRing(fSize,outerRadius,sigmaOuter);
    end
    myFilter = fspecial('gaussian',fSize,sigmaInner);
    myFilter((end+1)/2, (end+1)/2) = 0;% Set the center to 0
    myFilter = myFilter./sum(myFilter,'all'); %gaussian filter with 0 origin (but still sums to 1)
    

    if temporalFilter>0
        mov = imgaussfilt3(mov,[0,0,temporalFilter],"FilterSize",[1,1,temporalFilter*2+1]);
    end


    [h, w, frames] = size(mov);
    neighborAverage = zeros(h, w, frames);

    for i = 1:batchSize:frames
        idx = i:min(i+batchSize-1, frames);
        if sigmaOuter>0
            mov(:,:,idx) = mov(:,:,idx)-imfilter(gpuArray(mov(:,:,idx)),backgroundFilter,'circular'); %subtract background from everything
        end
        %after background subtract get neighbor average:
        neighborAverage(:,:,idx) = imfilter(gpuArray(mov(:,:,idx)),myFilter,"circular"); %circular so edges aren't brighter
    end



%     neighborAverage = imfilter(mov,myFilter,"circular"); %circular so edges aren't brighter
%     mov = mov-imfilter(mov,backgroundFilter,'circular');

    stdMain = std(mov,0,3);
    stdNeighbor = std(neighborAverage,0,3);
    
    mov = mov - mean(mov,3);
    neighborAverage = neighborAverage - mean(neighborAverage,3);
    correlationim = sum((mov.*neighborAverage),3)./(size(mov,3)-1);
    
    correlationim = correlationim./(stdMain.*stdNeighbor);
    correlationim = gather(correlationim);
end