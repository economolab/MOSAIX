function [fs,tmi,series] = load_time_series(params,DFOFFreq,regressBG,doCausal)
if ~exist("DFOFFreq",'var'); DFOFFreq=[]; end
if ~exist("regressBG","var"); regressBG = 0; end
if ~exist("doCausal","var"); doCausal = 1; end


tsFn = getFileName(params, 'TS');

f = load(tsFn);

for i = 1:numel(f.ts)
    f.ts{i} = f.ts{i}((params.vertSplit+1)*params.framesToDisregard+1:end, :);
    f.tsbg{i} = f.tsbg{i}((params.vertSplit+1)*params.framesToDisregard+1:end, :);
    f.tm{i} = f.tm{i}((params.vertSplit+1)*params.framesToDisregard+1:end);
end

% for i = 1:numel(ts)
%     ts{i} = ts{i}(params.framesToDisregard+1:end, :);
%     tsbg{i} = tsbg{i}(params.framesToDisregard+1:end, :);
%     tm{i} = tm{i}(params.framesToDisregard+1:end);
% end

% ts{1} = ts{1}(params.framesToDisregard+1:end, :);
% tsbg{1} = tsbg{1}(params.framesToDisregard+1:end, :);
% tm{1} = tm{1}(params.framesToDisregard+1:end);

fs = (1./mean(diff(f.tm{1})));
tmi = cell2mat(f.tm);
if size(tmi,1) ==1
    tmi = tmi';
end
tmi = tmi+1/fs;
series = cell2mat(f.ts);

if regressBG
    tsBG = cell2mat(f.tsbg);
    tsBG = tsBG-mean(tsBG);
    X = [ones(size(movingDisplacements, 1), 1), tsBG];

    beta = X \ series;

    % Calculate the fitted values for all neurons
    fittedValues = X * beta;

    % Compute the residuals for all neurons
    residuals = series - fittedValues;

    % Update the original data matrix with the residuals
    series = residuals;
end


if isempty(DFOFFreq)
    return
end
if (DFOFFreq>0)
    
    low = MySmoothFreqs(series,fs,DFOFFreq,doCausal);
    
    series = 100 * (series - low) ./ low; % Calculate the desired output
else
    for i = 1:nNeurons  % iterate over columns/neurons
        % Fit the decaying exponential model to the data
        
        f = fit((1:size(series,1))', medfilt1(series(:,i),5), 'exp1');  % 'exp1' models y = a*exp(b*x)
        
        % Evaluate the fitted model over the time vector
        fitted_exp = f.a .* exp(f.b .* (1:size(series,1))');
        
        series(:,i) = (series(:,i) - fitted_exp)*100./fitted_exp;
    end
end

end