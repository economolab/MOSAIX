function gaussRingFilt = getGaussianRing(fSize, radius,sigma)
    % Generate grid of coordinates with the center at the middle of the image
    [x, y] = meshgrid(1:fSize, 1:fSize);
    x_centered = x - (fSize+1)/2;
    y_centered = y - (fSize+1)/2;

    % Calculate distance from the radius of the ring
    distance = abs(sqrt(x_centered.^2 + y_centered.^2) - radius);

    % Generate Gaussian values
    gaussRingFilt = exp(-distance.^2 / (2*sigma^2));
    gaussRingFilt = gaussRingFilt./sum(gaussRingFilt,'all');
end
