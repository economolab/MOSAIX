function im = gaussFilt(im, c)

N = round(c*5);
if mod(N,2)==0
    N = N+1;
end
gaus2 = zeros(N,N);
midx = (N+1)/2;
midy = midx;

for x = 1:N
    for y = 1:N
        gaus2(x,y) = exp(-((x-midx)^2+(y-midy)^2)/(2*c^2));
    end
end
gaus2 = gaus2./sum(sum(gaus2));

for i = 1:size(im,3)
    temp = im(:,:,i);
    im(:,:,i) = 2.*conv2(temp, gaus2, 'same');
end
