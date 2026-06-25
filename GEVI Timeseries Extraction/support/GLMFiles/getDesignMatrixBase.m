function [designMatrixBase,correspondence,distancesToDrawnRois,spikes] = getDesignMatrixBase(param,neuronIndices)
    [fsStim, stim] = load_stim_series(param);
    tmStim = (1:length(stim)).*1/fsStim;
    
    [fsData,tm,data] = load_time_series(param);

    nCells = numel(neuronIndices);
    spikes = cell(1,nCells);
    for neuronNum = 1:nCells
        spikes{neuronNum} = getSpkIx(data(:,neuronNum), 3.5,0,param.highPassWindowSize);
    end
    stimOnset = [0;diff(stim)]>0;
    indicesOnset = floor(tmStim(stimOnset).*fsData);    
    indicesOnset = indicesOnset - param.framesToDisregard;
    indicesOnset(indicesOnset>length(tm)) = [];
    nStims = length(indicesOnset);

    [correspondence,distancesToDrawnRois] = load_stim_correspondence(param,1);
    distancesToDrawnRois = distancesToDrawnRois(neuronIndices,:);
    [nPatterns,nStimRois] = size(correspondence);
    nCells = length(spikes);
    designMatrixBase = zeros(length(tm),(nStimRois+nCells)); %(tm,stim)

    for j = 1:nStims
        correspondenceIndex = mod(j-1,nPatterns)+1;
        onsetIndex = indicesOnset(j);
        if onsetIndex<1
            continue
        end
        for stimRoiNum = find(correspondence(correspondenceIndex,:))
            designMatrixBase(onsetIndex,stimRoiNum) = 1;
        end
    end

    for neuronNum = 1:nCells
        for k = 1:length(spikes{neuronNum})
            designMatrixBase(spikes{neuronNum}(k),nStimRois+neuronNum) = 1;
        end
    end
end