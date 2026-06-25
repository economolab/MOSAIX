function data = getDataNaNStimResponse(data, fsData, tmOffset, stim, fsStim, stimCorrespondence, windowSize)
doPlot = 0;
if doPlot
    data2 = data;
end
% Nans data around each stim
% fsData is data sampling frequency
% fsStim is stim sampling frequency
% windowSize is number of samples on each side to NaN

% Find stim onset indices (rising edges)
onIndices = find(diff(stim) > 0);

% Number of neurons in the data
numNeurons = size(data, 2);

% Handle empty stimCorrespondence
if nargin < 5 || isempty(stimCorrespondence)
    % If stimCorrespondence is empty, assume all neurons correspond to all stims
    stimCorrespondence = true(1, numNeurons);
    numStimPatterns = 1;
else
    % Number of stim patterns
    numStimPatterns = size(stimCorrespondence, 1);
end

% Loop over each stim occurrence
for idx = 1:length(onIndices)
    % Determine the stim pattern index (cycles through available patterns)
    stimPatternIndex = mod(idx-1, numStimPatterns) + 1;
    
    % Get neurons associated with this stim pattern
    neurons = stimCorrespondence(stimPatternIndex, :)>0;
    
    % Stim index in stim data
    stimIndex = onIndices(idx);
    
    % Calculate the corresponding time of the stim
    t_stim = (stimIndex - 1) / fsStim;
    
    % Corresponding index in data
    dataIndex = round((t_stim-tmOffset) * fsData) + 1;
    
    % Define the range to set NaN around the stim
    startIndex = max(1, dataIndex - windowSize(1));
    endIndex = min(size(data, 1), dataIndex + windowSize(2));
    
    % Set data to NaN for the selected neurons in the specified window
    
    data(startIndex:endIndex, neurons) = NaN;
end
if doPlot
    nNeuron = size(data,2);
    sep = 8;
    figure;  hold on; 
    plot(data2+sep*(1:nNeuron),'k');plot(data+sep*(1:nNeuron));
end
end
