function ix = getSpkIx(dat, threshSD, do_plot,be_dumb,highPassWindowSize)
if (size(dat,1)==1)
    dat = dat';
end
if any(isnan(dat))
    dat = interpNan(dat);
end
if ~exist("windowSize","var") || isempty(highPassWindowSize)
    highPassWindowSize = 31;
end

lo = MySmooth(dat, highPassWindowSize,0);
dat = dat -lo;

pdat = [dat(dat>0); -dat(dat>0)];
mu = mean(pdat);
sd = std(pdat); %use std in positive direction to determine ~shot noise std
thresh = mu-threshSD*sd;

ix = find(dat<thresh & dat < circshift(dat,1) & dat < circshift(dat,-1));
% ix = find(dat(1:end-1)>thresh&dat(2:end)<thresh)+1;

if exist('be_dumb','var') && be_dumb
    return
end


Nbef = 2;
Naft = 3;
Npts = Nbef+Naft+1;


ix = ix(ix>Nbef);
ix = ix(ix<numel(dat)-Naft+1);
Npeaks = numel(ix);

snip = zeros(Npts, Npeaks);


for i = 1:Npeaks
    snip(:, i) = dat(ix(i)-Nbef:ix(i)+Naft);
end
meanwav = median(snip, 2);
meanwav = meanwav-mean(meanwav);
meanwav = meanwav./sum((abs(meanwav)));
a = xcorr(dat-mean(dat), -meanwav);
a = a(end-numel(dat)+1-Nbef:end-Nbef);
a = std(dat).*(a - mean(a))./std(a);

dat = a;

pdat = [dat(dat>0); -dat(dat>0)];
mu = mean(pdat);
sd = std(pdat);
thresh = mu-threshSD*sd;
ix = find(dat<thresh & dat < circshift(dat,1) & dat < circshift(dat,-1));

if ~exist('do_plot', 'var') || ~do_plot
    return
end
fig = figure(23);
clf(fig, 'reset');
ax1 = subplot(1, 5, 1, 'Parent',fig); 
ax2 = subplot(1, 5, 2:5,'Parent',fig);
plot(ax1,snip); hold(ax1,'on'); plot(ax1,mean(snip, 2), 'k', 'LineWidth', 2); 
plot(ax2,dat); hold(ax2,'on'); plot(ax2,a);
