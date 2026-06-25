function Pipeline_getTimeSeries(params)

ExpInfoFn = getFileName(params, 'ExpInfo');
RoiFn = getFileName(params, 'ROI');

load(ExpInfoFn,'fn');
if exist(RoiFn,'file')
    load(RoiFn,'roi');
else
    roi = [];
end
fn.pth = regexprep(params.pth, '(\\|/)', filesep);
fn.si = params.imtoken;

[ts, tm, tsbg] = getAllTimeSeries(fn, roi, params);


wsfn = getFileName(params, 'wsfn');

if exist(wsfn,'file')==2 && isscalar(tm) && tm{1}(1)==0
    offset = getTMOffset(wsfn);
end


% tm = tm(params.framesToDisregard:end, :, :); 
% ts = ts(params.framesToDisregard:end, :, :);
% tsbg = tsbg(params.framesToDisregard:end, :, :);

didReg = params.doRegistration;

BGRemovalRadius = params.BGRemovalRadius;

save(getFileName(params, 'TS'), 'ts', 'tm', 'tsbg', 'didReg', 'BGRemovalRadius', '-v7.3');

end

function offset = getTMOffset(wsfn)
    offset = 0; return;
    data = ws.loadDataFile(wsfn);
    index = find(contains(lower(data.header.AIChannelNames),'all rows'));
    if isempty(index)
        warning('could not identify which channel contained all rows frame times');
        offset = 0;
        return
    end
    
    allFramesExposed = data.sweep_0001.analogScans(:,index)>1;
    allFramesExposed(find([0; diff(allFramesExposed)>0],1,'last'):find([0; diff(allFramesExposed)<0],1,'last')) = 0; %0 out last one since it isn't real
    offset = find(allFramesExposed,1)/data.header.AcquisitionSampleRate; %offset in seconds
    sample_rate = data.header.StimulationSampleRate;
end