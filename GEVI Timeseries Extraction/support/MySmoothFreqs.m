function out = MySmoothFreqs(x, fs, fc, causal)
if ~exist('causal','var'); causal=1; end

if size(x,1)==1 && size(x,2)>1
    x=x';
end

% Calculate N based on fs and fc
N = round(0.68 * fs / fc);

if N<=0
    out = x;
    return;
end

if mod(N, 2)==0
    N = N+1;
end

kern = gausswin(N);
if causal
    kern(1:floor(N/2)) = 0; %causal
end

kern = kern./sum(kern);

Ncol = size(x, 2);
Nel = size(x, 1);

out = zeros(Nel, Ncol);
for j = 1:Ncol
    out(:, j) = conv(x(:, j), kern, 'same');
    for i = 1:ceil(N/2)    
        out(i, j) = mean(x(1:i, j));
    end
end
