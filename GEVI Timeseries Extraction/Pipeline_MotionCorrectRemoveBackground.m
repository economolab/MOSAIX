function Pipeline_MotionCorrectRemoveBackground(params)
addpath("..\NoRMCorre")

ExpInfoFn = getFileName(params, 'ExpInfo');
regMovFn = getFileName(params, 'RegBGMov');
regInfoFn = getFileName(params, 'RegInfo');
load(ExpInfoFn,'fn');
fn.pth = regexprep(params.pth, '(\\|/)', filesep);
fn.si = params.imtoken;
fnIn = fullfile(fn.pth,fn.si{1});
batchSize = 200; startFrameTemplate = 800; medianFiltRange = 1;

fr = FrameReader(fnIn);

if startFrameTemplate>fr.maxFrames
    startFrameTemplate = round(fr.maxFrames/2);
end

if exist(regMovFn,'file')
    if exist(regInfoFn,'file')
        regInfoPrev = load(regInfoFn);
        if params.BGRemovalRadius == regInfoPrev.bgRadius && length(imfinfo(regMovFn)) == fr.maxFrames
            return
        else
            delete(regMovFn);
        end
    else
        delete(regMovFn);
    end
end


frames = fr.getFramesByIndexes(max(startFrameTemplate,1),min(startFrameTemplate+batchSize-1,fr.maxFrames));

ogType = class(frames);

options_rigid = NoRMCorreSetParms('d1',size(frames,1),'d2',size(frames,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',batchSize,'iter',0,'print_msg',0,'usegpu',1);
options_rigid.iter = 1;
options_rigid.correct_bidir = 0;

[~,~,template] = normcorre(frames,options_rigid);
 
fr.pointer=1;
clear futureObj
shifts = [];
tstart = tic;
f = waitbar(0,"Running registration + BG subtraction");
fTIF = Fast_BigTiff_Write(regMovFn,1,0);
if medianFiltRange > 1
    numExcludeStart = ceil(medianFiltRange/2)-1;
    framesForMedFilt = repmat(fr.getFramesByIndexes(1,1),[1,1,numExcludeStart]);
else
    framesForMedFilt = [];
end
while (~fr.isDone)
    [M_final,shiftsNew] = normcorre(fr.getFrames(batchSize),options_rigid,template); %runs on GPU as rolling ball runs on cpu
    shifts = [shifts;shiftsNew]; %#ok<AGROW> 
    if exist("futureObj",'var')
        tempSave(fetchOutputs(futureObj),ogType,fTIF);
        futureObj = [];
    end

    if params.BGRemovalRadius>0
        futureObj = parfeval(backgroundPool,@rollingBall,1,M_final,params.BGRemovalRadius,medianFiltRange,framesForMedFilt);
        framesForMedFilt = M_final(:,:,end-medianFiltRange+2:end);
    else
        tempSave(M_final,ogType,fTIF);
    end


    time = toc(tstart);
    waitbar(fr.pointer/fr.maxFrames,f,sprintf('Reg + BG. %d / %d + %d seconds',round(time),round(time),round(time./fr.pointer*(fr.maxFrames-fr.pointer))));
end

if exist("futureObj",'var')
    tempSave(fetchOutputs(futureObj),ogType,fTIF);
end

if params.BGRemovalRadius>0 && ~isempty(framesForMedFilt)
    lastFrames = rollingBall(framesForMedFilt(:,:,ceil(medianFiltRange/2):end),params.BGRemovalRadius,medianFiltRange,framesForMedFilt);
    tempSave(lastFrames,ogType,fTIF);
end

fTIF.close()

shifts = [shifts.shifts];

delete(f);

if size(shifts,1)~=size(shifts,3) || size(shifts,1)~=1
    keyboard; %what are 1st and third dimensions
else
    shifts = squeeze(shifts);
end

bgRadius = params.BGRemovalRadius;
if exist(regInfoFn,'file')
    save(regInfoFn,"shifts","bgRadius",'-append');
else
    save(regInfoFn,"shifts","bgRadius");
end

%===============
function tempSave(data,ogType,fTIF)

data = cast(data,ogType);
for i =1:size(data,3)
    fTIF.WriteIMG(data(:,:,i));
end