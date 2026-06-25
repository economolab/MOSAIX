function [sample_rate, stim] = load_stim_series(param,stdThresh,askForThresh)

    if (~exist("stdThresh","var") || isempty(stdThresh)); stdThresh = []; end
    if (~exist("askForThresh","var") || isempty(askForThresh)); askForThresh = 1; end


    fullfn = getFileName(param, 'OptoData');
    
    if exist(fullfn,"file")
        f = load(fullfn,"sample_rate","stim");
        sample_rate = f.sample_rate;
        stim = f.stim;
        return
    end

    wsfn = getFileName(param, 'wsfn');

    if exist(wsfn,'file')==2
        [sample_rate,stim] = loadWsfnData(wsfn);
    elseif (stdThresh<=0)
        sample_rate = [];
        stim = [];
        return
    else
        stim = detectLightArtifact(param,stdThresh,askForThresh);
        [sample_rate,~,~] = load_time_series(param); 
    end

    if (size(stim,1)==1)
        stim = stim';
    end
    if (size(sample_rate,1)==1)
        sample_rate = sample_rate';
    end
    save(fullfn,"stim","sample_rate");
end

function [sample_rate,stim] = loadWsfnData(wsfn)
    data = ws.loadDataFile(wsfn);
    stim = data.sweep_0001.analogScans(:,3)>1;
    if ~contains(lower(data.header.AIChannelNames(3)),'stim')
        warning("potentially wrong index for AIChannel for wavesurfer data")
    end
    if (stim(end))
        stim(end+1)=0;
    end
    sample_rate = data.header.StimulationSampleRate;
end

function [sample_rate, stim] = load_stim_series_helper(param)

    matFile = dir(fullfile(param.pth, '*.mat'));
    a=load(fullfile(param.pth,matFile(1).name));
    stim = a.analogInData(:,1);
    sample_rate = a.analogInSamplingRate;
    return

    % f = load(getFileName(param, 'WSData'), 'wsdata');
    % stim = f.wsdata.scan{1}.digitalScans>1;
    % sample_rate = f.wsdata.meta.AcquisitionSampleRate;
end