function myDotProdIm = dotProductWithBackground(param)

if isfield(param,'framesPerBatch')
    framesPerBatch = param.framesPerBatch;
else
    framesPerBatch = 400;
end

tsFn = getFileName(param, 'TS');
load(tsFn, 'tsbg');

if isfield(param,"offset")
    offset = param.offset;
else
    offset = [];
end

if param.doRegistration
    Pipeline_MotionCorrectRemoveBackground(param)
end

fn.pth = regexprep(param.pth, '(\\|/)', filesep);
fn.si = param.imtoken;


if numel(fn.si)>1
    throw('error-- not supported for multiple files');
end

if param.doRegistration
    imReader = FrameReader(getFileName(param,'RegBGMov'));
else
    imReader = FrameReader(fullfile(fn.pth,fn.si{1}));
end

MovavgFn = getFileName(param, 'Movavg');
load(MovavgFn,'meanim');

myDotProdIm = zeros(size(meanim));
meanim = rollingBall(meanim,param.BGRemovalRadius);
minVal = prctile(meanim,5,"all");
if minVal ==0
    minVal = 5; %idk arbitrary
end
minValMask = meanim<minVal;
meanim = double(meanim);
meanimDenom = meanim; 
meanimDenom(minValMask) = inf;
while ~imReader.isDone(framesPerBatch,true)
        frames = imReader.getFrames(framesPerBatch);
    
        if isempty(offset) && prctile(frames,1,'all')>95
            offset = 100;
        elseif isempty(offset)
            offset = 0;
        end
        frames = double(frames);
        frames = frames-offset;
        dffPixels = (frames-meanim)./meanimDenom;
        dffPixels = medfilt3(dffPixels,[1,1,5]);
        dffPixels = medfilt3(dffPixels,[5,5,1]);
        runningAvOfAll = squeeze(mean(frames,[1,2]));
        runningAvOfAll = medfilt1(runningAvOfAll,5);
        dffAll = (runningAvOfAll-median(runningAvOfAll))./median(runningAvOfAll);
        myDotProdIm = myDotProdIm + sum(dffPixels.*reshape(dffAll,1,1,[]),3);
end

end