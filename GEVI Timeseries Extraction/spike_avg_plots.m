function spike_avg_plots(params, threshSD,numIndividualTrialsToPlot,neuronIndices,distanceMat,includeBootstrap,stimCorrespondence)

doFigureVersion = 1;

if ~exist("threshSD",'var');threshSD = 3.5;end
if ~exist("includeBootstrap",'var');includeBootstrap = 0;end
if ~exist("numIndividualTrialsToPlot",'var') || isempty(numIndividualTrialsToPlot)
    numIndividualTrialsToPlot = 0;
end


fig = figure();
clf(fig);
ax = gca; hold(ax,'on');
hpFreq = 1;
[fsData,tm,data] = load_time_series(params,hpFreq);
[fsStim, stim] = load_stim_series(params,0);
nextStimTime = ceil(median(diff(find(diff(stim)>0)))/fsStim*fsData); %num frames til next stim occurs
nextStimTime =nextStimTime -2; %subtract two frames (I guess?) for delay to spike


% data = data(1:2000,:);
if exist('neuronIndices','var')
    data = data(:,neuronIndices);
    stimCorrespondence = stimCorrespondence(:,neuronIndices);
    distanceMat = distanceMat(neuronIndices,neuronIndices);
end
nNeurons = size(data, 2);


fontSize = 40;

% data = data - MySmooth(medfilt1(data,150),25);



dataNoDirectOpto = getDataNaNStimResponse(data,fsData,tm(1),stim,fsStim,stimCorrespondence,[5,15]);

data = interpNan(data);

num_neurons = size(data,2);

xPadding = 15; % Define the x-axis offset
numBefore = 10;
numAfter = 30;
window_size = numBefore+1+numAfter; % Define the size of the window around spike times
xOffset = xPadding + window_size;

if doFigureVersion
    smoothingFactor = 15;
else
    smoothingFactor=1;
end

dataNoSpikes = data;
spxs = {}; idxs = {};
for i = 1:num_neurons
    spxs{i} = getSpkIx(data(:,i), threshSD, 0, params.highPassWindowSize);
    idxs{i} = spxs{i} + (-numBefore:numAfter); idxs{i} = idxs{i}';
    idxs{i}(idxs{i}<1 | idxs{i}>size(data,1)) = [];
    dataNoSpikes(idxs{i},i) = nan;
end



numSpikesAllCells = zeros(num_neurons,1);
avg_data = getAvgData(data);
avg_data_no_direct_opto = getAvgData(dataNoDirectOpto);

if includeBootstrap
    numSamplesForDistribution = 10000;percentiles =[2.5,97.5];
    bootstrapFN = getFileName(params,'BootstrapInfo');
    bootstrapDat = [];
    if exist(bootstrapFN,'file')==2
        bootstrapDat = load(bootstrapFN);
        if numel(bootstrapDat.numSpikesAllCells)~=numel(numSpikesAllCells) || ~all([bootstrapDat.smoothingFactor==smoothingFactor,bootstrapDat.numSamplesForDistribution==numSamplesForDistribution,all(bootstrapDat.numSpikesAllCells==numSpikesAllCells),bootstrapDat.window_size==window_size,all(bootstrapDat.percentiles==percentiles)])
            bootstrapDat = [];
        end
    end
    if ~isempty(bootstrapDat)
        bootstrapPercentiles = bootstrapDat.bootstrapPercentiles;
    else
        f = msgbox(sprintf('Running bootstrap (takes ~2 min for 39 neurons Θ(n^2))'));
        bootstrapPercentiles = getBootstrapPercentiles(MySmooth(data,smoothingFactor,1),numSamplesForDistribution,numSpikesAllCells,window_size,percentiles);
        save(bootstrapFN,"bootstrapPercentiles","numSpikesAllCells","numSamplesForDistribution","window_size","percentiles","smoothingFactor");
        delete(f);
    end
end


% Sigmoid function for color blending
% sigmoid = @(x) 1 ./ (1 + exp((-x+30)/8));
sigmoid = @(x) 1 ./ (1 + exp((-x+50)/16));
totalIndividualsToPlot = sum(cellfun(@(indCell) size(indCell,3), avg_data(:,2)));
if totalIndividualsToPlot>500
    warning('plotting %d individual trials',totalIndividualsToPlot)
    % keyboard;
end

if numIndividualTrialsToPlot>0
    sep = max(cat(3,avg_data{:,2}),[],'all');
else
    maxedData = cat(3,avg_data{:,1});
    for j = 1:size(maxedData,2)
        maxedData(:,j,j) = 0;
    end
    ranges = max(maxedData)-min(maxedData);
    sep = prctile(ranges(:),97);
    sep = sep*1.1;
    clear maxedData;
end

addToPlot(avg_data)

xlim([0,xOffset*num_neurons]);
ylims =[-sep/2,sep*num_neurons];
ylim(ylims);

nextStimRectangleWidth = window_size-nextStimTime-numBefore;

yOffsets = (0:num_neurons-1)*sep;
lineAlpha = 0.2;
for j = 1:num_neurons
    currentXOffset = (j - 1) * xOffset;
    currentYOffset = (j - 1) * sep;
    xline(currentXOffset+numBefore+1,'Alpha',lineAlpha)
    if includeBootstrap
        for i = 1:num_neurons
            if i==j
                continue
            end
            x = [currentXOffset, currentXOffset + window_size];
            yVals = squeeze(yOffsets(i) + bootstrapPercentiles(:, i, j));
            
            
            % Create x and y vectors for the polygon to fill
            xPolygon = [x, fliplr(x)];  % Repeat x values in reverse for the return path
            yPolygon = [yVals(1), yVals(1), yVals(2), yVals(2)];
            
            % Use fill to create a shaded region
            fill(ax, xPolygon, yPolygon, 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        end
        % yline(currentYOffset + bootstrapPercentiles(:,j),'Alpha',lineAlpha);
    end
    yline(currentYOffset,'Alpha',lineAlpha);
    rectangle(Position=[currentXOffset,currentYOffset-sep/4,window_size,sep])
    if ~doFigureVersion && nextStimRectangleWidth>0
        rectangle(Position=[currentXOffset+nextStimTime+numBefore,ylims(1),nextStimRectangleWidth,ylims(2)-ylims(1)],FaceColor='r',FaceAlpha=0.03,EdgeColor="none");
    end

end

if doFigureVersion
    scalebarPos1 =[-5,0,2,2]; %#ok<UNRCH>
    rectangle(ax,Position=scalebarPos1,FaceColor=[0 0 0 0],EdgeColor='none',FaceAlpha=1)
    text(ax,scalebarPos1(1),scalebarPos1(2)+scalebarPos1(4)/2,sprintf('%d%%',scalebarPos1(4)),HorizontalAlignment='right',FontSize=50)
    ms = 100;
    scalebarPos2 =[0,-1,ms*fsData/1000,0.1];
    rectangle(ax,Position=scalebarPos2,FaceColor=[0 0 0 0],EdgeColor='none',FaceAlpha=1)
    text(ax,scalebarPos2(1)+scalebarPos2(3)/2,scalebarPos2(2),sprintf('%d ms',ms),HorizontalAlignment='center',VerticalAlignment='top',FontSize=50)
    axis(ax,'off');
end
xlim(ax,[-5,xOffset*(nNeurons)]);

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

set(fig,'color','w');

setappdata(fig,"axes",ax)
setappdata(fig,"fs",fsData)
setappdata(fig,"window_size",window_size)
setappdata(fig,"xPadding",xPadding)

% fig.SizeChangedFcn = @figResize;
% set(zoom(ax),'ActionPostCallback',@(x,y) figResize(fig));
% set(pan(fig),'ActionPostCallback',@(x,y) figResize(fig));
% addlistener(ax, 'XLim', 'PostSet', @(src, ev) figResize(fig));
    function avg_data = getAvgData(dataIn)
        avg_data = cell(num_neurons,2);
        for ii = 1:num_neurons
            spike_idx = getSpkIx(data(:,ii), threshSD); %get spikes based on data with direct opto
            numSpikesAllCells(ii) = length(spike_idx);
            temp_data = zeros(window_size, num_neurons, numSpikesAllCells(ii));
            
            for jj = 1:numSpikesAllCells(ii)
                start_idx = max(1, spike_idx(jj) - numBefore);
                end_idx = min(size(dataIn, 1), spike_idx(jj) + numAfter);
        
        
                %in case the max/min take effect
                temp_data_start_idx = 1+start_idx-(spike_idx(jj) - numBefore);
                temp_data_end_idx = window_size-((spike_idx(jj) + numAfter)-end_idx);
                temp_data(temp_data_start_idx:temp_data_end_idx, :, jj) = dataIn(start_idx:end_idx, :);
            end
            % temp_data = temp_data-mean(temp_data,1);
            
            if numIndividualTrialsToPlot && numSpikesAllCells(ii)>0
                numSpikesAllCells(ii) = size(temp_data,3); 
                indices = 1:numSpikesAllCells(ii);
                if numSpikesAllCells(ii) > numIndividualTrialsToPlot
                    indices = round(linspace(1, size(temp_data,3), numIndividualTrialsToPlot));
                end
                avg_data{ii,2} = temp_data(:,:,indices);
            end
            % avg_data{ii,1} = squeeze(mean(temp_data, 3,'omitnan'));
            avg_data{ii,1} = squeeze(median(temp_data, 3,'omitnan'));
        end
    end

    function addToPlot(avg_data)

        for ii = 1:num_neurons
            currentXOffset = (ii - 1) * xOffset;
            for jj = 1:num_neurons
                % Calculate the current plot's offsetst
                currentYOffset = (jj - 1) * sep;
                color = [1-sigmoid(distanceMat(ii,jj)), 0, 0]; % Red to Black transition
                spikeAverage = squeeze(-avg_data{ii,1}(:,jj));
                spikeAverage = spikeAverage - mean(spikeAverage,'omitnan');
                lineWidth = 1;
                factor = 1;
                if ii==jj
                    factor = (1/max(spikeAverage(:)).*sep/2);
                elseif any(spikeAverage(numBefore+2:end)>0.2)
                    lineWidth = 2;
                end
                if ii~=jj
                    spikeAverage = MySmooth(spikeAverage,smoothingFactor,1);
                end
        
                if numIndividualTrialsToPlot
                    data = squeeze(-avg_data{ii,2}(:, jj, :));
                    plot(ax,(1:window_size) + currentXOffset, data.*factor + currentYOffset, 'DisplayName', sprintf('Neuron %d,%d', ii, jj));
                end
                
                plot(ax,(1:window_size) + currentXOffset, spikeAverage.*factor + currentYOffset, 'Color',color, 'LineWidth',lineWidth, 'DisplayName', sprintf('Neuron %d,%d', ii, jj));
            end
        end
    end
end

function customColormap = getColormap(sigmoid,xRange,distanceMat)
    nColors = 256; % Number of colors in the colormap
    % Generate colormap
    customColormap = zeros(nColors, 3); % Initialize colormap matrix
    xRange = linspace(xRange(1), xRange(2), nColors); % Range of x values for sigmoid function
    for i = 1:nColors
        blendRatio = sigmoid(xRange(i));
        customColormap(i, :) = [1-blendRatio, 0, 0]; % Red to black
    end

    if ~exist('distanceMat','var') || isempty(distanceMat)
        return
    end

    figure;imagesc(distanceMat);
    % Apply the custom colormap
    colormap(customColormap);

    % Assuming you have a matrix 'data' to visualize
    colorbar; % Add a colorbar to visualize the colormap

    
end


function figResize(src, ~)
end





% [~,~,data] = load_time_series(params);
% % data = data(1:2000,:);
% if exist('neuronIndices','var')
%     data = data(:,neuronIndices);
%     distanceMat = distanceMat(neuronIndices,neuronIndices);
% end
% data = detrend(data, 3,'omitnan');
% data = interpNan(data);
% 
% num_neurons = size(data,2);
% 
% xPadding = 2; % Define the x-axis offset
% window_size = 50; % Define the size of the window around spike times
% xOffset = xPadding + window_size;
% 
% avg_data = cell(num_neurons,2);
% 
% for i = 1:num_neurons
%     spike_idx = getSpkIx(data(:,i), threshSD);
%     num_spikes = length(spike_idx);
%     temp_data = zeros(window_size, num_neurons, num_spikes);
% 
%     for j = 1:num_spikes
%         start_idx = max(1, spike_idx(j) - floor(window_size/2));
%         end_idx = min(size(data, 1), spike_idx(j) + floor(window_size/2)-1);
% 
% 
%         %in case the max/min take effect
%         temp_data_start_idx = 1+start_idx-(spike_idx(j) - floor(window_size/2));
%         temp_data_end_idx = window_size-((spike_idx(j) + floor(window_size/2)-1)-end_idx);
%         temp_data(temp_data_start_idx:temp_data_end_idx, :, j) = data(start_idx:end_idx, :);
%     end
%     temp_data = temp_data-mean(temp_data,1);
% 
%     if numIndividualTrialsToPlot && num_spikes>0
%         numSpikes = size(temp_data,3); 
%         neuronIndices = 1:numSpikes;
%         if numSpikes > numIndividualTrialsToPlot
%             neuronIndices = round(linspace(1, size(temp_data,3), numIndividualTrialsToPlot));
%         end
%         avg_data{i,2} = temp_data(:,:,neuronIndices);
%     end
%     avg_data{i,1} = squeeze(mean(temp_data, 3));
% end
% 
% % Sigmoid function for color blending
% sigmoid = @(x) 1 ./ (1 + exp((-x+30)/8));
% 
% for i = 1:num_neurons
%     for j = 1:num_neurons
%         % Calculate the current plot's offsetst
%         currentYOffset = (j - 1) * sep;
%         currentXOffset = (i - 1) * xOffset;
%         if any(isnan(avg_data{i,1}))
%             continue
%         end
% 
%         if numIndividualTrialsToPlot
%             data = squeeze(-avg_data{i,2}(:, j, :));
%             plot(ax,(1:window_size) + currentXOffset, data + currentYOffset, 'DisplayName', sprintf('Neuron %d,%d', i, j));
%         end
% 
%         dist = distanceMat(i,j);
%         blendRatio = sigmoid(dist); % Adjust this factor to control the transition sharpness
%         color = [1-blendRatio, 0, 0]; % Red to Black transition
%         plot(ax,(1:window_size) + currentXOffset, squeeze(-avg_data{i,1}(:,j)) + currentYOffset, 'Color',color, 'LineWidth',3, 'DisplayName', sprintf('Neuron %d,%d', i, j));
%     end
% end    
% hold off; % Release the plot
% colorLims = [0,100];
% colormap(ax,getColormap(sigmoid,colorLims));
% colorbar(ax); clim(ax,colorLims);


% old way of plotting in tiled layout (I don't like because you can't zoom)
% t = tiledlayout(fig, num_neurons, num_neurons);
% t.TileSpacing = 'compact';
% t.Padding = 'tight';
% 
% % Initialize a variable to store axes for linking
% allAxes = gobjects(num_neurons^2, 1);
% axIdx = 1; % Index for storing axes handles
% 
% for i = 1:num_neurons
%     for j = 1:num_neurons
%         ax = nexttile(t, (num_neurons - i) * num_neurons + j);
%         plot(ax, -avg_data(:, i, j));
%         
%         % Store the axis handle
%         allAxes(axIdx) = ax;
%         axIdx = axIdx + 1;
% 
%         % Remove axis labels for inner tiles
%         if i ~= 1
%             ax.XTickLabel = [];
%         end
%         if j ~= 1
%             ax.YTickLabel = [];
%         end
%     end
% end
% 
% % Link all axes for consistent zooming and panning
% linkaxes(allAxes, 'xy');
% for i = 1:num_neurons
%     vec = true(1,num_neurons);
%     vec(i)=0;
%     ylim(allAxes(1+0*i),[min(-avg_data(:, i, vec),[],"all"),max(-avg_data(:, i, vec>0),[],"all")])
% end


