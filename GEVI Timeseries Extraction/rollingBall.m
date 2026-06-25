function image = rollingBall(image, radius, medianFiltRange, framesFromOtherBatch)
    if exist('framesFromOtherBatch','var') && ~isempty(framesFromOtherBatch)
        image = cat(3,framesFromOtherBatch,image);
    end

    if ~exist('medianFiltRange','var')
        medianFiltRange = 1;
    end

    if (radius<=10) %copied from Imagej
        shrinkfactor = 1;
    elseif (radius<=30)
        shrinkfactor = 2;
    elseif (radius<=100)
        shrinkfactor = 4;
    else
        shrinkfactor = 8;
    end
    
    bg = imresize(image,1/shrinkfactor,"bilinear");
    if radius<=2 %idk it's messed up for rad==2
        se = offsetstrel('ball',round(radius/shrinkfactor),round(radius/shrinkfactor),0); % figure;imagesc(se.Offset)
        bg = vectorizedGpuOpening(bg,se.Offset);
    else
        se = offsetstrel('ball',round(radius/shrinkfactor),round(radius/shrinkfactor)); % figure;imagesc(se.Offset)
        bg = imopen(bg,se);
    end

    rows = size(image,1); cols = size(image,2);
    bg = imresize(bg,[rows,cols],'bilinear');
    if medianFiltRange>1
        numExcludeStart = ceil(medianFiltRange/2)-1;
        numExcludeEnd = floor(medianFiltRange/2);
        ogType = class(bg);
        bg = medfilt1(single(bg),medianFiltRange,[],3);
        
        image = image(:,:,numExcludeStart+1:end-numExcludeEnd);
        bg = bg(:,:,numExcludeStart+1:end-numExcludeEnd);
        
        bg = cast(bg,ogType);
    end

    image = image-bg;
end

% This is faster when "n-Number of nonflat line-shaped structuring elements
% used to approximate the shape" argument of imopen is 0
function outputImage = vectorizedGpuOpening(inputImage, structElem)
    batchSize = 50; %idk faster doing it on smaller batches which I wouldn't have guessed
    % Pad input image to handle border effects
    [SErows, SEcols] = size(structElem);
    padSizeRow = floor(SErows / 2);
    padSizeCol = floor(SEcols / 2);
    paddedImage = padarray(inputImage, [padSizeRow, padSizeCol], 'replicate', 'both');
    
    paddedImage = single(paddedImage);
    [ix,iy] = find(isfinite(structElem));
    structElemVals = structElem(isfinite(structElem));
    ix = ix-padSizeRow; iy = iy - padSizeCol;
    % Preallocate an array to hold all shifted images
    [pRows, pCols,pImages] = size(paddedImage,[1,2,3]);
    resultImage = inf(pRows, pCols, pImages,'single'); % Use nnz for number of non-zero elements in structElem
    outputImage = inputImage;
    for batchStart = 1:batchSize:size(inputImage,3)
        batchEnd = min(batchStart+batchSize-1,size(inputImage,3));
        resultImageBatch = gpuArray(resultImage(:,:,batchStart:batchEnd));
        paddedImageBatch = gpuArray(paddedImage(:,:,batchStart:batchEnd));
        ix = gpuArray(ix); iy = gpuArray(iy);
        for i = 1:numel(ix)
            % Create shifted image and store in the 3D matrix
            resultImageBatch = min(resultImageBatch,circshift(paddedImageBatch, [ix(i), iy(i)])-structElemVals(i));
        end
        paddedImageBatch = resultImageBatch;
        for i = 1:numel(ix)
            resultImageBatch = max(resultImageBatch,circshift(paddedImageBatch, [ix(i), iy(i)])+structElemVals(i));
        end
        outputImage(:,:,batchStart:batchEnd) = gather(resultImageBatch(padSizeRow+1:end-padSizeRow, padSizeCol+1:end-padSizeCol,:,:));
    end
end
