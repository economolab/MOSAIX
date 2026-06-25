function out = MySmooth(x, N, causal, dim)
if ~exist('causal','var'); causal=1; end

if size(x,1)==1 && size(x,2)>1 && ~exist('dim','var')
    x=x';
end

if ~exist('dim','var'); dim=1; end

if N<=0
    out = x;
    return;
end

Ncol = size(x, 2);
Nel = size(x, 1);

if mod(N, 2)==0
    N = N+1;
end

kern = gausswin(N);
if causal
    kern(1:floor(N/2)) = 0; %causal
end

kern = kern./sum(kern);



out = convn(x,kern,'same');
% for j = 1:Ncol
%     out(:, j) = conv(x(:, j), kern, 'same');
%     for i = 1:ceil(N/2)    
%         out(i, j) = mean(x(1:i, j));
%     end
% end


