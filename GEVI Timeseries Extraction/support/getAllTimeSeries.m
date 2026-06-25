function [ts, tm, tsbg]= getAllTimeSeries(fn, mask, params)
    mask = logical(mask);
    
    if ~isfield(fn,'si')
        fn.si={""}; %#ok<STRSCALR> 
    end
    
    if isfield(params,"offset")
        offset = params.offset;
    else
        offset = [];
    end


    Nfiles = numel(fn.si);
    
    ts = cell(Nfiles, 1);
    tm = cell(Nfiles, 1);
    tsbg = cell(Nfiles, 1);

    if params.doRegistration
        Pipeline_MotionCorrectRemoveBackground(params)
    end
    
    for i = 1:numel(fn.si)
        noRegReader = FrameReader(fullfile(fn.pth,fn.si{i}));
        if params.doRegistration
            if i>1
                throw('error-- not supported for multiple files')
            end
            regReader = FrameReader(getFileName(params,'RegBGMov'));
            if sum(mask,'all')>0
                [ts{i}, ~] = getTimeSeries(regReader, mask,offset,params);
            end
            
            if exist(getFileName(params,'TS'),'file')
                contents = load(getFileName(params,'TS'));
                tsbg{i} = contents.tsbg{i}; tm{i} = contents.tm{i};
                continue;
            end
            [~,tsbg{i}] = getTimeSeries(noRegReader, mask,offset,params);
        else
            [ts{i}, tsbg{i}] = getTimeSeries(noRegReader, mask,offset,params);
        end

        tm{i} = getTm(noRegReader);
    end
end

function [ts, tsbg] = getTimeSeries(imReader, mask,offset,params)
    if ~exist("offset","var")
        offset = [];
    end
    Nframes = imReader.maxFrames;
    Nroi = size(mask,3);
    framesPerBatch=400;
    ts= zeros(Nframes,Nroi);
    tsbg= zeros(Nframes,1);
    neuropilRange = [];
    if exist('params','var') && isfield(params,'neuropilRange'); neuropilRange = params.neuropilRange; end
    
    while ~imReader.isDone(framesPerBatch,true)
        start = imReader.pointer;
        frames = imReader.getFrames(framesPerBatch);
        finish = imReader.pointer-1;
    
        if isempty(offset) && prctile(frames,1,'all')>95
            offset = 100;
        elseif isempty(offset)
            offset = 0;
        end

        frames = frames-offset;
        
    
        for k =1:Nroi
            ts(start:finish,k) = getTimeSeriesBinaryMask(frames,mask(:,:,k),neuropilRange,0,sum(mask,3));
        end
        if numel(mask)>0
            tsbg(start:finish) = getTimeSeriesBinaryMask(frames,sum(mask,3)==0,neuropilRange,0);
        else
            tsbg(start:finish) = squeeze(mean(frames,[1,2]));
        end
    end
end

function tmi = getTm(imreader)
%need to do tsbg without BG subtraction to be able to detect stims
    if ~isempty(imreader.metadata)
        tmi = [imreader.metadata.tm];
        tmi = tmi/1000000;
        tmi = tmi - min(tmi);
        while 1/mean(diff(tmi))<5
            tmi = tmi/1000;
        end
        return;
    end

    fnFullNoReg  = imreader.fileID;
    if isfolder(fnFullNoReg)
        folder = fnFullNoReg;
    else
        folder = fileparts(fnFullNoReg); % Get the directory part of the filename
    end
    csvFiles = dir(fullfile(folder, '*.csv'));
    if length(csvFiles)==1
        frameTimes = table2array(importfile_CSV_frametimes(fullfile(folder,csvFiles.name)));
        frameTimes=frameTimes-frameTimes(1);
        tmi=frameTimes/1000000; %convert to seconds
        return
    end


    try
        imageDescriptionCellArray = {imfinfo(fn).ImageDescription};
        tmi = cellfun(@(in) getMetaVal(in,'meta.header.timeBof'),imageDescriptionCellArray);
        tmi = tmi/1000000000000; %ps to s
        tmi = tmi - min(tmi);
    catch
        tmi = 1:Nframes;
    end
end