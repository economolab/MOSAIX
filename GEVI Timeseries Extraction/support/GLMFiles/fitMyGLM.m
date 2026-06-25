function [B,pVals] = fitMyGLM(opts,designMatrix,neuronTimeSeries,nBasisFunctions,neuronIdxs,dmdStimWeightsThisCell)
    if opts.HowToUseOwnSpikes == 0 %do nothing
        designMatrix(:,neuronIdxs) = 0; %don't include this neuron's spikes in the fit
    elseif opts.HowToUseOwnSpikes == 1 %include as beta coeffs
    elseif opts.HowToUseOwnSpikes == 2 %exclude that time
        spikes = getSpkIx(neuronTimeSeries, 3.5, 0, opts.param.highPassWindowSize);
        spikeWindowFrames = reshape((spikes+(-opts.spikeWindowToRemove(1):opts.spikeWindowToRemove(2))),[],1);
        spikeWindowFrames(spikeWindowFrames<1 | spikeWindowFrames>size(designMatrix,1))=[];
        designMatrix(spikeWindowFrames,:) = [];
        neuronTimeSeries(spikeWindowFrames) = [];
    end
    
    if opts.weightDMDStimByDistance
        designMatrix = designMatrix.*dmdStimWeightsThisCell;
    end
    B = lasso(designMatrix, neuronTimeSeries,'Alpha',opts.Alpha,'Lambda',opts.Lambda);
    
    if opts.DoStats
        pVals =  nonparametricBootstrap(opts,designMatrix,neuronTimeSeries,opts.Lambda,nBasisFunctions);
    else
        pVals=0;
    end



    % function [isSignificant,pVals] = myWaldTest()
    %     mdl = fitglm(designMatrixIndividual, neuronTimeSeries);
    %     goodVec = any(mdl.CoefficientCovariance~=0);
    %     nBeta = size(designMatrixIndividual,2)+1;
    %     r = cell(nCells+nStimRois,1);
    %     covMats = cell(nCells+nStimRois,1);
    %     mdlEst =mdl.Coefficients.Estimate(goodVec);
    %     R = cell(1,nCells+nStimRois);
    %     for i = 1:nCells+nStimRois
    %         idxVec = false(1,nBeta);
    %         idxStart = nBasisFunctions*(i-1)+1;
    %         idxs = idxStart:idxStart+nBasisFunctions-1;
    %         idxVec(idxs)=1;
    %         idxVec = idxVec(goodVec);
    %         idxs = find(idxVec);
    %         r{i} = mdlEst(idxVec)';
    %         R{i} = zeros(sum(idxVec),sum(goodVec));
    %         for j = 1:numel(idxs)
    %             R{i}(j,idxs(j)) = 1;
    %         end
    %         covMats{i} = mdl.CoefficientCovariance(goodVec,goodVec);
    %     end
    %     [isSignificant,pVals] = waldtest(r,R,covMats);
    % end
end


function pValues = nonparametricBootstrap(opts,designMat,response,lambda,nBasisFunctions)
    % Number of bootstrap samples
    nBoot = 500;
    n = size(designMat, 1);
    bootstrapCoefs = zeros(size(designMat,2)/nBasisFunctions, nBoot);
    idx = zeros(n,nBoot);
    for i = 1:nBoot
        idx(:,i) = randsample(n, n, true);
    end
    tic
    parfor i = 1:nBoot
        % Resample with replacement
        yb = response(idx(:,i));
        designMatb = designMat(idx(:,i),:);

        % Fit the Lasso GLM model to the bootstrap sample
        [B_boot, ~] = lasso(designMatb, yb,'Alpha',opts.Alpha, 'Lambda', lambda);
        bootstrapCoefs(:, i) = mean(reshape(B_boot,nBasisFunctions,[]));
    end
    toc
    prctiles = zeros(2,size(bootstrapCoefs,1));
    prctiles(1,:) = sum(bootstrapCoefs>=0,2)/size(bootstrapCoefs,2);
    prctiles(2,:) = sum(bootstrapCoefs<=0,2)/size(bootstrapCoefs,2);
    pValues = min(prctiles)*2;
end
