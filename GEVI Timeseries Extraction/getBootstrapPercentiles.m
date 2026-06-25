function bootstrapPercentiles = getBootstrapPercentiles(data, numSamplesForDistribution, numSamplesToAverage,windowSize, percentiles)
    % bootstrapPercentiles(:,j,i) is percentiles of neuron j meaned over numSamplesToAverage(i)
    % Number of time points
    numTimePoints = size(data, 1);
    
    % Number of neurons
    numNeurons = size(data, 2);
    bootstrapPercentiles = zeros(numel(percentiles),numNeurons,length(numSamplesToAverage));
    if isscalar(numSamplesToAverage)
        resampleIndices = randi(numTimePoints-windowSize, [numSamplesToAverage, numSamplesForDistribution]);
    else
        resampleIndices = randi(numTimePoints-windowSize, [max(numSamplesToAverage), numSamplesForDistribution]);
    end
    resampleIndices = permute(resampleIndices,[1,3,2]) + (1:windowSize);
    for j = 1:numNeurons
        neuronData = data(:,j);
        sampledData = neuronData(resampleIndices);
        if isscalar(numSamplesToAverage)
            dat = squeeze(mean(sampledData));
            dat = squeeze(max(dat)); %max over each window
            bootstrapPercentiles(:,j) = prctile(dat, percentiles);
        else
            for i = 1:numNeurons
                dat = squeeze(mean(sampledData(1:numSamplesToAverage(i),:,:))); %mean however many windows we had 
                dat = dat-mean(dat);
                bootstrapPercentiles(percentiles<50,j,i) = prctile(squeeze(min(dat)), percentiles(percentiles<50));
                bootstrapPercentiles(percentiles>50,j,i) = prctile(squeeze(max(dat)), percentiles(percentiles>50));
            end
        end
    end
    
    % Calculate the percentiles for each time point
end