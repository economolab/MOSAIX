function plotGLMValidation(betaCoef,designMatrix,data,tm,stimWeights,spikes,spikeWindowToRemove)
    figure;ax = gca(); hold(ax,"on");
    ySepPredicted = 10;
    colors = brewermap(3,'Dark2');
    dataFit = zeros(size(data));
    for neuronNum = 1:size(betaCoef,3)
        spikeWindowFrames = reshape((spikes{neuronNum}+(-spikeWindowToRemove(1):spikeWindowToRemove(2))),[],1);
        spikeWindowFrames(spikeWindowFrames<1 | spikeWindowFrames>size(data,1))=[];

        yCoord = (neuronNum-1)*ySepPredicted;
        betas = betaCoef(:,:,neuronNum);
        stimWeightsThisCell = reshape(repmat(stimWeights(:,neuronNum)',[size(designMatrix,2)/size(stimWeights,1),1]),1,[]);
        designMatrixIndividual = designMatrix.*stimWeightsThisCell;
        
        dataFit(:,neuronNum) = glmval(betas(:),designMatrixIndividual,'identity','constant','off');
        dataFit(spikeWindowFrames,neuronNum) = nan;
        a = plot(ax,tm,yCoord-data(:,neuronNum),DisplayName='raw data',Color=colors(1,:));
        b = plot(ax,tm,yCoord-dataFit(:,neuronNum),DisplayName='fit data',Color=colors(2,:)) ;
        c = plot(ax,tm,yCoord+dataFit(:,neuronNum)-data(:,neuronNum),DisplayName='residuals',Color=colors(3,:)); 
        yline(yCoord,LineWidth=0.2,Alpha=0.4)
    end
    legend(ax,[a,b,c]); %just so it only has 3 entries
end
