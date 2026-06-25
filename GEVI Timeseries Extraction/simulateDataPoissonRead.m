function noised_data = simulateDataPoissonRead(meanVal,gain,readNoise,area,sampleLength)
    % Generate a dataset following a Poisson distribution
    
    poisson_data = poissrnd(meanVal.*gain, sampleLength, area)/gain;% Per pixel measurement
    readNoise = normrnd(0, readNoise, sampleLength, area); %Per pixel measurement
    
    poisson_data = sum(poisson_data,2)./size(poisson_data,2); %ROI measurement
    readNoise = sum(readNoise,2); %ROI measurement

    noised_data = poisson_data + readNoise;
end