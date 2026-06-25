function [stimsOut,fsStim] = getThisNeuronsStims(param,roiIndex,stimCorrespondence,stim)    
    if ~exist('stimCorrespondence','var') || isempty(stimCorrespondence)
        roiFile = load(getFileName(param, 'ROI'));
        stimCorrespondence = roiFile.stimCorrespondence;
    end
    if ~exist('stim','var') || isempty(stim)
        [fsStim,stim] = load_stim_series(param,0);
    else
        fsStim = [];
    end
    numStims = numel(find(diff(stim)>0));
    numRepeatsNeeded = ceil(numStims/size(stimCorrespondence,1));
    stimCorrespondence = repmat(stimCorrespondence,[numRepeatsNeeded,1]);
    stimCorrespondence = stimCorrespondence(1:numStims,:);

    roiStims = find(stimCorrespondence(:,roiIndex));

    stimsOut = false(size(stim));
    timesOn = find(diff([0;stim])>0);
    timesOff = find(diff([0;stim])<0);
    for k = 1:length(roiStims)
        idxStart = timesOn(roiStims(k));
        idxEnd = timesOff(find(timesOff>idxStart,1));
        stimsOut(idxStart:idxEnd) = 1;
    end
end



% roiIndex = j;
% dmdCorrespondence = roiFile.stimCorrespondence;