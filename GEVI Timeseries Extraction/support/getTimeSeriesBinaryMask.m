function time_series = getTimeSeriesBinaryMask(im, mask, neuropilRange,normalize,allROIMask,timeFilterBGSize)
    if ~exist('neuropilRange',"var") || isempty(neuropilRange)
        neuropilRange=[0 0];
    end
    if ~exist('timeFilterBGSize','var')
        timeFilterBGSize = 1;
    end
    % Compute the mean for each frame
    im = reshape(im, [], size(im,3));
    time_series = mean(im(mask > 0.5, :), 1);

    if neuropilRange(2)-neuropilRange(1)>0
        seBG = strel('disk', neuropilRange(2));
        seIgnore = strel('disk', neuropilRange(1));
        maskOuter = ~imerode(double(~mask),seBG)-~imerode(double(~mask),seIgnore);
        if exist("allROIMask",'var') && ~isempty(allROIMask) && ... 
                sum(maskOuter&(~allROIMask),'all')/sum(maskOuter,'all')>=0.4
            %only use non-roi for background (unless there's less than 40%
            %non-roi in that region in which case just use the roi too)
            maskOuter=maskOuter&(~allROIMask);
        end
    %     mask = reshape(mask, height, width);
        bg = MySmooth(mean(im(maskOuter > 0.5, :), 1),timeFilterBGSize);
        if any(size(bg)~=size(time_series))
            bg = bg';
        end
        time_series = time_series - bg;
    end
    
    %default to normalize so only if passed and <= do we return
    if ~exist("normalize","var") || normalize<=0
        return
    end
    % Compute the median and normalize the time_series
    mu = median(time_series);
    time_series = 100 * (time_series - mu) ./ mu;
end
