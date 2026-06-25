function light_artifact_flag = detectLightArtifact(param, std_threshold,askForThresh,excludePercent)
if ~exist("std_threshold",'var') || isempty(std_threshold);std_threshold = 6;end
if ~exist("askForThresh",'var') || isempty(askForThresh);askForThresh = 0;end
if ~exist("excludePercent",'var');excludePercent = 10;end

tsFn = getFileName(param, 'TS');
if ~exist(tsFn,'file')
    Pipeline_getTimeSeries(param)
end
load(tsFn, 'tsbg');
if isscalar(tsbg) && iscell(tsbg)
    tsbg = tsbg{1};
else
    error('unexpected format')
end
if size(tsbg,2) ==1
    tsbg = tsbg';
end

if all(tsbg==0)
    tsbg = getFrameMeans(param);
end

diffFrames = [0, diff(tsbg)];
toCalculateSTDOn = diffFrames<prctile(diffFrames,100-excludePercent) & diffFrames>prctile(diffFrames,excludePercent);
standardDev = std(diffFrames(toCalculateSTDOn));
diffFrames = diffFrames./standardDev; 

if askForThresh
    std_threshold = getThreshold(diffFrames);
end

light_artifact_flag_on = diffFrames>std_threshold;
light_artifact_flag_on(light_artifact_flag_on&circshift(light_artifact_flag_on,-1))=0;
light_artifact_flag_off = diffFrames<-std_threshold;
light_artifact_flag_off(light_artifact_flag_off&circshift(light_artifact_flag_off,1))=0;

light_artifact_flag = cumsum(light_artifact_flag_on)-cumsum(light_artifact_flag_off);
light_artifact_flag = light_artifact_flag-median(light_artifact_flag)>0;

if askForThresh
    return
end

timeBetweenStims = sort(diff(find(light_artifact_flag_on)));
if (sum(light_artifact_flag)<10 || length(timeBetweenStims)<10 || std(timeBetweenStims(3:end-1))>1)
    std_threshold = getThreshold(diffFrames);
    light_artifact_flag_on = diffFrames>std_threshold;
    light_artifact_flag_on(light_artifact_flag_on&circshift(light_artifact_flag_on,1))=0;
    light_artifact_flag_off = diffFrames<-std_threshold;
    light_artifact_flag_off(light_artifact_flag_off&circshift(light_artifact_flag_off,1))=0;
    
    light_artifact_flag = cumsum(light_artifact_flag_on)-cumsum(light_artifact_flag_off);
    light_artifact_flag = light_artifact_flag-median(light_artifact_flag)>0;
end

if false %code to correct
    std_threshold = 1;
    figure;plot(light_artifact_flag_on); hold on; plot(light_artifact_flag_off); yyaxis right; plot(tsbg)
    figure;plot(diff(find(light_artifact_flag_on)));
    figure;plot(diff(find(light_artifact_flag_off)));
    stdDiffs = std(diff(find(light_artifact_flag_on)));
    medianDiffs = median(diff(find(light_artifact_flag_on)));
    while true
        indiciesOn = find(light_artifact_flag_on);
        badIndicies = [0, diff(indiciesOn)<medianDiffs-stdDiffs*4];
        if all(~badIndicies)
            break
        end
        light_artifact_flag_on(indiciesOn(find(badIndicies,1)))=0;
    end
    while true
        indiciesOn = find(light_artifact_flag_off);
        badIndicies = [0, diff(indiciesOn)<medianDiffs-stdDiffs*4];
        if all(~badIndicies)
            break
        end
        light_artifact_flag_off(indiciesOn(find(badIndicies,1)))=0;
    end
    light_artifact_flag = cumsum(light_artifact_flag_on)-cumsum(light_artifact_flag_off);
    light_artifact_flag = light_artifact_flag-median(light_artifact_flag)>0;
end

end

function yout = getThreshold(data)
    f = figure;
    ax = gca();
    plot(ax, data);
    title(ax, 'Click inside the plot area to activate selection mode.');
    
    f.ButtonDownFcn = @(src, event) axesClicked(ax);
    ax.ButtonDownFcn = @(src, event) axesClicked(ax);
    y_select = 0;

    uiwait(f);
    yout = y_select;

    function axesClicked(src)
        % Disable further axis click detections to prevent re-entry
        src.ButtonDownFcn = '';
        
        % Inform the user to select a point
        title(src, 'Now select a point. Click once to choose.');
        
        % Wait for user to select a point
        [~, y_select] = ginput(1);
        
        delete(f);
    end
end



% thresh2 = standardDev*(std_threshold*1.5);
% light_artifact_flag_on = diffFrames>thresh1  | (diffFrames-circshift(diffFrames,-1)>thresh2) | (diffFrames+circshift(diffFrames,-1)>thresh2);
% figure;
% title("Stim detected (note as long as after the first ignored frames its ok)");
% ax1 = subplot(3, 1, 1); plot(diffFrames); hold on; yline(thresh1); plot(diffFrames>thresh1)
% ax2 = subplot(3, 1, 2); plot(diffFrames-circshift(diffFrames,-1)); hold on; yline(thresh2); plot(diffFrames-circshift(diffFrames,-1)>thresh2)
% ax3 = subplot(3,1,3); plot(diffFrames+circshift(diffFrames,-1)); hold on; yline(thresh2); plot(diffFrames+circshift(diffFrames,-1)>thresh2)
% linkaxes([ax1,ax2,ax3],'x')
% keyboard  %for editing if not working as best as possible

function frameMeans = getFrameMeans(param)
ExpInfoFn = getFileName(param, 'ExpInfo');

tempmask.mask = [];
meanIm = load(getFileName(param,"Movavg"),'meanim').meanim;
selectROI2(meanIm, @drawrectangle,@setMask,"Select region for detecting stim");
uiwait(); %have to wait until figure done--- it's annoying but return vars of this function are dependent
stimDetectionMask = tempmask.mask;

stimDetectionMask1D = stimDetectionMask(:);
if sum(stimDetectionMask1D,"all")==0
    stimDetectionMask1D(:)=1;
end

load(ExpInfoFn,'fn');
fn.pth = regexprep(param.pth, '(\\|/)', filesep);
fn.si = param.imtoken;
frameMeans = [];
for i = 1:numel(fn.si)
    imReader = FrameReader(fullfile(fn.pth,fn.si{i}));
    
    Nframes = imReader.maxFrames;
    framesPerBatch=400;
    frameMeansFile= zeros(1,Nframes);
    
    while ~imReader.isDone(framesPerBatch,true)
        start = imReader.pointer;
        frames = imReader.getFrames(framesPerBatch);
        finish = imReader.pointer-1;
        
        % Element-wise multiplication with mask
        frames = reshape(frames, [], size(frames, 3));
        frameMeansFile(start:finish) = mean(frames(stimDetectionMask1D,:),1);
    end
    frameMeans(end+1:end+Nframes) = frameMeansFile;
end

function setMask(maskOut,~)
    tempmask.mask = maskOut; % Assign the mask variable when the app is completed
    uiresume();
end
end


