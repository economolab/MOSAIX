function simple_plot(param,filt,threshSD,dispThresh,sep,axesToPlot,flipDFF,indicesToPlot,dmdROICorrespondence,stopBandFreq,regressOutMotion,colors,hpFreq)
    % figure;axesToPlot(1) = gca();
    if exist('axesToPlot','var') && ~isempty(axesToPlot)
        axesToPlotTS =axesToPlot(1); axesToPlotMotion = axesToPlot(2);
        cla(axesToPlotTS); cla(axesToPlotMotion);
    end

    indicesToPlot= sort(indicesToPlot);
    if ~exist("filt",'var')
        filt = 1;
    end
    if ~exist("threshSD",'var')
        threshSD = 3;
    end


    movingWorking =0;
    regInfoFn = getFileName(param,'RegInfo');
    if exist(regInfoFn,'file')
        regInfo = load(regInfoFn);
        movingDisplacements = regInfo.shifts;
        movingDisplacements = movingDisplacements(param.framesToDisregard+1:end,:);
    else
        movingWorking =0;
    end

    if ~exist("axesToPlotTS",'var') || isempty(axesToPlotTS)
        fig=figure(333);
        clf(fig);
        axesToPlotTS = axes(fig);

        if movingWorking
            ax1 = subplot(5,1,1:4,'Parent',axesToPlotTS); % Create subplot for the main plot
            ax2 = subplot(5,1,5,'Parent',axesToPlot); % Create subplot for the moving displacement plot
        end
    else
        ax1 = axesToPlotTS;
        if movingWorking
            axesToPlotTS.Layout.Row(2) = max(axesToPlotMotion.Layout.Row,[],'all')-1;
            axesToPlotMotion.Visible='on';
            ax2 = axesToPlotMotion;
        else
            axesToPlotTS.Layout.Row(2) = max(axesToPlotMotion.Layout.Row,[],'all');
            axesToPlotMotion.Visible='off';
        end
    end
    hold(ax1,'on');
    if ~exist("flipDFF",'var')
        flipDFF = -1;
    end
    flipMult = (flipDFF>0)*-2+1;
% 
    [fsData,tm,data] = load_time_series(param,hpFreq);
    data  = data(:,indicesToPlot);
    num_neurons = length(indicesToPlot);


    if isempty(tm) || (max(tm,[],"all")-min(tm,[],"all")==0)
        tm = (0:length(data)-1)*1/1000;
        title(axesToPlot,"NO METADATA FOR TIMESTAMPS FOUND")
    end
%     data = detrend(data, 3,'omitnan');
    
    
    if movingWorking && exist("regressOutMotion",'var') && ~isempty(regressOutMotion) && regressOutMotion>0
        % Loop through each neuron
        X = [ones(size(movingDisplacements, 1), 1), movingDisplacements];
        if regressOutMotion>1
            X = [X, movingDisplacements.*movingDisplacements, movingDisplacements(:,1).*movingDisplacements(:,2)]; %#ok<AGROW> 
        end

        beta = X \ data;
    
        % Calculate the fitted values for all neurons
        fittedValues = X(:,2:end) * beta(2:end,:);
    
        % Compute the residuals for all neurons
        residuals = data - fittedValues;
    
        % Update the original data matrix with the residuals
        data = residuals;
    end

    if exist('filterFreq','var') && ~isempty(stopBandFreq) && stopBandFreq>0
        BW = 2; % Bandwidth of the notch filter (in Hz)
        Wn = [(stopBandFreq-BW)/(fsData/2) (stopBandFreq+BW)/(fsData/2)];
        [b, a] = butter(4, Wn, 'stop'); % Design the Butterworth filter
        
        for i = 1:size(data,2)
            data(:,i) = filtfilt(b, a, data(:,i)); % Apply the Butterworth filter
        end
        title(axesToPlot,"APPLYING BUTTERWORTH TO REMOVE " + string(stopBandFreq) + " FREQUENCY")
    end

    isNanIdx = find(isnan(data));
    data(isNanIdx) = data(isNanIdx-1);


    for i=1:num_neurons
        % tm = (61:length(tm)+60).*1/fsData;
        plot(ax1,tm,sep*(i-0.5)+flipMult*MySmooth(data(:,i),filt),'Color',colors(indicesToPlot(i),:));
        
        highPassWindowSize = 1/10*fsData;

        spike_idx = getSpkIx(data(:,i), threshSD,0,highPassWindowSize);
        if ~isempty(spike_idx)
            plot(ax1,tm(spike_idx),sep*(i),'.k');
        end
    end

    ylabel(ax1,'-\Delta F/F')
    grid(ax1,'on'); yticks(ax1,0:10:num_neurons*sep);

    labels = string([-sep/2,repmat(-sep/2+10:10:sep/2,1,num_neurons+1)]);
    labels(strcmp(labels,string(sep/2)))="(-)"+string(sep/2);
    yticklabels(ax1,labels)
    
    xlim(ax1,[0,max(tm(:))]);

    if movingWorking
        plot(ax2,tm,movingDisplacements);
        hold(ax2,'on');
        yline(ax2,dispThresh)
        ylabel(ax2,'Moving Displacement')    
        ylim(ax2,'auto');
        zoom(ax2,'reset');

        % linkaxes([ax1,ax2],'x');
    end

    [fsStim,stim] = load_stim_series(param);


    
    if sep>0
        ylim(ax1,[-sep/2,sep*num_neurons+sep/2]);
    else
        ylim(ax1,'auto');
    end
    ylims = ylim(ax1);
    zoom(ax1,'reset');
    if isempty(dmdROICorrespondence)
        area(ax1,(0:length(stim)-1).*1/fsStim,(stim>0.2)*ylims(2), 'FaceColor',[0 0 1], 'FaceAlpha',.3,'EdgeAlpha',0);
        return
    end
    stimPatterns = dmdROICorrespondence(:,indicesToPlot);
    if ~isempty(stim)

        stimPatternsIndex=1;
        onIndices = find(diff(stim)>0)'; offIndicies = find(diff(stim)<0)';
        width = mean(offIndicies-onIndices)/fsStim;
        x = (mean([offIndicies;onIndices])-width/2)/fsStim;
        for idx = 1:length(onIndices)
            for j = find(stimPatterns(stimPatternsIndex,:))
                rectangle(ax1,'Position',[x(idx),(j-1)*sep,width,sep],'FaceColor',[0 0 1 0],'EdgeColor',[0,0,0,0.2],'FaceAlpha',0.2)
            end
            stimPatternsIndex=mod(stimPatternsIndex,size(stimPatterns,1))+1;
        end
    end


end


