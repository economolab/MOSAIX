function GLMPlot(opts)
% GLMPlot Plots Generalized Linear Model (GLM) for neural data analysis
%
% Arguments:
% param - Structure containing parameters for data and analysis
% DistanceMat - Matrix of distances between neurons (optional)
% NeuronIndices - Indices of neurons to include in the analysis (optional)
% BasisFunctions - Basis functions for GLM (default: depends on UseExpBasis)
% UseExpBasis - whether to use exponential basis when no basis functions supplied
% HPFreq - High-pass filter frequency (default: 0.5 Hz)
% RegressOwnSpikes - Boolean to regress own spikes (default: 0)
% Lambda - Regularization parameter for LASSO. -1 is calculate with cross validation (default: 0)
% DoStats - Boolean to perform bootstrap statistical analysis (default: true if lambda is 0)
%
% Example usage:
% GLMPlot(param, 1:15, distanceMat, eye(20), 0.5, 0, 0, true)
%
% Notes:
% - If basisFunctions is not provided or is empty, it defaults to an identity matrix.
% - If HPFreq is not provided or is empty, it defaults to 0.5 Hz.
% - If regressOwnSpikes is not provided or is empty, it defaults to false (0).
% - If lambda is not provided or is empty, it defaults to 0 (no regularization).
% - If doStats is not provided or is empty, it defaults to true if lambda is 0, otherwise false.
%
% The function processes the provided neural data and fits a GLM to analyze
% the relationships between neuron activity and various predictors.    

    arguments
        opts.param
        opts.DistanceMat double = -1;
        opts.NeuronIndices double = [];
        opts.BasisFunctions double = [];
        opts.nBasisFunctions double = 50;
        opts.HPFreq double = 0.5
        opts.HowToUseOwnSpikes double = 2;  %0-- do nothing 1-- include with beta 2-- exclude time around
        opts.spikeWindowToRemove (1,2) double = [3,6]; % if HowToUseOwnSpikes ==2 excludes this many frames before/after 
        opts.Lambda double = 0.00575
        opts.crossValidationFolds double = 20;
        opts.DoStats double = 1;
        opts.Alpha double = 1;
        opts.weightDMDStimByDistance = 1; %whether to make it so further DMD distances are less likely to be used
        opts.numToShiftBasisFunctions double = 3; %how many frames to circshift design matrix to have things before spikes
        opts.plotOwnSpikes double = 0;
    end
    %% get opts
    if isempty(opts.BasisFunctions)
        opts.BasisFunctions = eye(opts.nBasisFunctions);
    else
        opts.BasisFunctions = opts.BasisFunctions;
    end
    nBasisFunctions = size(opts.BasisFunctions,2);
    
    [fsData,tm,data] = load_time_series(opts.param,opts.HPFreq);
    if isempty(opts.NeuronIndices); opts.NeuronIndices=1:size(data,2); end
    data = data(:,opts.NeuronIndices);
    nCells = size(data,2);
    if isempty(opts.DistanceMat)
        opts.DistanceMat = zeros(nCells);
    else
        opts.DistanceMat = opts.DistanceMat(opts.NeuronIndices,opts.NeuronIndices);
    end
    
    %% check if done this calculation before
    saveFn = getFileName(opts.param, 'GLM');
    if exist(saveFn,'file')
        glmFileContents = load(saveFn);
        if opts.Lambda == -1
            resetLam=true;
            opts.Lambda = glmFileContents.opts.Lambda;
        else
            resetLam=false;
        end
        if isequaln(glmFileContents.opts,opts)
            load(saveFn,'betaCoef','designMatrix','tm','stimWeights','spikes','opts','responses','fsData','plotBold','distancesToDrawnRois');
            plotGLMValidation(betaCoef,designMatrix,data,tm,stimWeights,spikes,opts.spikeWindowToRemove);
            plotResponses(betaCoef,responses,fsData,plotBold,opts.DistanceMat,distancesToDrawnRois,opts.numToShiftBasisFunctions,opts.plotOwnSpikes)
            return;
        end
        if resetLam
            opts.Lambda=-1;
        end
    end
    %%
    
    nTimePoints = size(opts.BasisFunctions,1);
        

    [designMatrixBase,correspondence,distancesToDrawnRois,spikes] = getDesignMatrixBase(opts.param,opts.NeuronIndices);
    
    [nPatterns,nStimRois] = size(correspondence);
    
    getDMDStimIdx = @(roiNum) (1 + (roiNum-1)*(nBasisFunctions)):(1 + (roiNum-1)*(nBasisFunctions)) + nBasisFunctions-1;
    getNeuronIdxs = @(roiNum) (1 + (roiNum-1 + nStimRois)*(nBasisFunctions):1 + (roiNum-1 + nStimRois)*(nBasisFunctions)+nBasisFunctions-1);
    
    designMatrix = getDesignMatrixFromBasisFunctions(designMatrixBase,opts.BasisFunctions,opts.numToShiftBasisFunctions);
    designMatrix(length(tm)+1:end,:) = [];    

    betaCoef = zeros(nBasisFunctions,nStimRois+nCells,nCells);
    stimWeights = ones(nStimRois+nCells,nCells);
    responses = zeros(nTimePoints,nStimRois+nCells,nCells);
    plotBold = false(nStimRois+nCells,nCells);
    pValues = zeros(nStimRois+nCells,nCells);

    if opts.weightDMDStimByDistance
        sigmoidDMDDistance = @(x) 1 ./ (1 + exp((x-25)/6)); %figure; plot(sigmoidDMDDistance(0:100)); %#ok<UNRCH> 
        for neuronNum = 1:nCells
            %arbitrary sigmoid function
            dmdStimWeightsThisCell = sigmoidDMDDistance(distancesToDrawnRois(neuronNum,:));
            dmdStimWeightsThisCell(dmdStimWeightsThisCell<0.03) = 0; %exclude things small
            stimWeights(1:nStimRois,neuronNum) = dmdStimWeightsThisCell;
        end
    end
    if opts.Lambda==-1
        
        designMatrixIndividual = designMatrix;
        if opts.HowToUseOwnSpikes == 0 %do nothing
            designMatrixIndividual(:,neuronIdxs) = 0; %don't include this neuron's spikes in the fit
        elseif opts.HowToUseOwnSpikes == 1 %include as beta coeffs
        elseif opts.HowToUseOwnSpikes == 2 %exclude that time
            spikes = getSpkIx(neuronTimeSeries, 3.5, 0, opts.param.highPassWindowSize);
            spikeWindowFrames = reshape((spikes+(-opts.spikeWindowToRemove(1):opts.spikeWindowToRemove(2))),[],1);
            spikeWindowFrames(spikeWindowFrames<1 | spikeWindowFrames>size(designMatrixIndividual,1))=[];
            designMatrixIndividual(spikeWindowFrames,:) = [];
            neuronTimeSeries(spikeWindowFrames) = [];
        end
        
        if opts.weightDMDStimByDistance
            dmdStimWeightsThisCell = reshape(repmat(stimWeights(:,neuronNum)',[nBasisFunctions,1]),1,[]);
            designMatrixIndividual = designMatrixIndividual.*dmdStimWeightsThisCell;
        end
        [B,FitInfo] = lasso(designMatrixIndividual, neuronTimeSeries,'Alpha',opts.Alpha, 'NumLambda',50,'CV',opts.crossValidationFolds,Options=statset(UseParallel=true));
        B = B(:,FitInfo.IndexMinDeviance);
        opts.Lambda = FitInfo.LambdaMinDeviance;
        
    end

    for neuronNum = 1:nCells
        
        dmdStimWeightsThisCell = reshape(repmat(stimWeights(:,neuronNum)',[nBasisFunctions,1]),1,[]);
        [B,pValues(:,neuronNum)] = fitMyGLM(opts,designMatrix,data(:,neuronNum),nBasisFunctions,feval(getNeuronIdxs, neuronNum),dmdStimWeightsThisCell); %#ok<FVAL>
        betaCoef(:,:,neuronNum) = reshape(B,nBasisFunctions,[]);
        
    end
    for neuronNum = 1:nCells
        for stimIdx = 1:size(betaCoef,2)
            responseVals = squeeze(opts.BasisFunctions*betaCoef(:,stimIdx,neuronNum)*stimWeights(stimIdx,neuronNum));
            responses(:,stimIdx,neuronNum) = responseVals;
        end
    end
    plotBold = pValues<0.05;
    save(saveFn,'betaCoef','designMatrix','tm','stimWeights','spikes','responses','fsData','plotBold','distancesToDrawnRois','opts');
    plotGLMValidation(betaCoef,designMatrix,data,tm,stimWeights,spikes,opts.spikeWindowToRemove);
    plotResponses(betaCoef,responses,fsData,plotBold,opts.DistanceMat,distancesToDrawnRois,opts.numToShiftBasisFunctions,opts.plotOwnSpikes)
    return

    

    function plotDataAlignedToSameStim()
        designMatrixBaseStims = designMatrixBase(:,1:nStimRois);
        ySepTemp = 0;
        numToCheck = 5;
        dexToPlot = 1;
        firstWithSomething = find(sum(designMatrixBaseStims,2),numToCheck);
        shifts = cell(numToCheck,1);
        for i = 1:numToCheck
            dex = firstWithSomething(i);
            shifts{i} = find(designMatrixBaseStims*designMatrixBaseStims(dex,:)'==sum(designMatrixBaseStims(dex,:)))-dex;
        end
        temp = cellfun(@size,shifts,'UniformOutput',false);temp = cat(1,temp{:});
        temp = temp(:,1);
        badIndices =temp ~= mode(temp);
        shifts(badIndices) = [];
        shifts = round(mean([shifts{:}]'));
        %shifts(end+1) = size(data,1);
        
        figure; ax = gca(); hold(ax,'on');
        numToPlot = length(shifts)-1;
        
        colors = brewermap(numToPlot,'RdBu');
        colors = colors./min(1,sum(colors')'); % make the really pale colors less so
        

        dataFit = zeros(size(data));
        dataLength = floor(mean(diff(shifts)));
        sameData = zeros(dataLength,numToPlot);
        for i = 1:numToPlot
            x = shifts(i)+1:shifts(i)+dataLength;
            sameData(:,i) = data(x,dexToPlot);
            plot(ax,-sameData(:,i)+ySepTemp*i,Color=colors(i,:))
            % betas = betaCoef(:,:,neuronIdx);
            % dataFit(:,neuronIdx) = glmval(betas(:),designMatrix,'identity','constant','off');
            % plot(ax,dataFit(x,dexToPlot)-data(x,dexToPlot)+ySepTemp*i,Color=colors(i,:))
        end
    end
end





% %**********unfinished
% function plotBold = CalcPlotBoldResponses(responses)
%     plotBold = zeros(size(responses));
% end
% 
% function plotBold = CalcPlotBoldPValues(pVals,nTimePoints)
%     nBootstrap = 50000;
%     nToMult = nTimePoints+1;
%     randVals = rand(nToMult,nBootstrap);
%     thresh = prctile(prod(randVals),5);
%     totalComparisons = size(pVals,2)*size(pVals,3);
%     plotBold = prod(pVals)<thresh/totalComparisons; %bonferoni
% end

