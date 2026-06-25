function plotResponses(betaCoef,responses,fsData,plotBold,distanceMatNeurons,distancesToDrawnRois,numToShiftBasisFunctions,includeOwnSpikes)
    fontSize = 40;
    additionalLineWidth = 1.5;
    ySep = 2;
    xSep = 3;
    nTimePoints = size(responses,1);
    stimSpikeSep = nTimePoints;

    tmLags = (0:nTimePoints-1)/fsData;
    figure; ax = gca();
    hold(ax,"on");

    % Sigmoid function for color blending
    sigmoid = @(x) 1 ./ (1 + exp((-x+50)/16));
    nStimRois = size(distancesToDrawnRois,2);
    nCells = size(distancesToDrawnRois,1);

    for stimIdx = 1:size(betaCoef,2)
        xCoord = (stimIdx-1)*(nTimePoints+1+xSep);
        xCoord = xCoord + stimSpikeSep*(stimIdx>nStimRois);
        xCoord = xCoord/fsData;
        for neuronNum = 1:nCells
            if stimIdx<=nStimRois
                dist = distancesToDrawnRois(neuronNum,stimIdx);
            else
                dist = distanceMatNeurons(stimIdx-nStimRois,neuronNum);
                if stimIdx-nStimRois==neuronNum && ~includeOwnSpikes
                    continue
                end
            end
            blendRatio = sigmoid(dist);
            color = [1-blendRatio, 0, 0]; % Red to Black transition
            yCoord = (neuronNum-1)*ySep;
            lineWidth = 0.5 + plotBold(stimIdx,neuronNum)*additionalLineWidth;
            plot(ax,xCoord+tmLags,yCoord-responses(:,stimIdx,neuronNum),LineWidth=lineWidth,Color=color);
        end
    end
    for neuronNum = 1:size(betaCoef,3)
        yCoord = (neuronNum-1)*ySep;
        yline(yCoord,LineWidth=0.2)
    end

    for stimIdx = 1:size(betaCoef,2)
        xCoord = (stimIdx-1)*(nTimePoints+1+xSep);
        xCoord = xCoord + stimSpikeSep*(stimIdx>nStimRois);
        xCoord = xCoord/fsData;
        xline(xCoord + tmLags(numToShiftBasisFunctions+1),LineWidth=0.2,Alpha=0.2);
    end

    colorRange = [0,100];
    nColors = 256; % Number of colors in the colormap
    % Generate colormap
    customColormap = zeros(nColors, 3); % Initialize colormap matrix
    colorLims = linspace(colorRange(1), colorRange(2), nColors); % Range of x values for sigmoid function
    for i = 1:nColors
        blendRatio = sigmoid(colorLims(i));
        customColormap(i, :) = [1-blendRatio, 0, 0]; % Red to black
    end
    colormap(ax,customColormap);
    clim(ax,colorRange);
    c = colorbar;
    c.Ruler.TickLabelFormat='%g um';
    c.Ruler.FontSize = floor(fontSize/2);
end