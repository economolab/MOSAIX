function output = normalize_local_contrast(image, blockRadius, meanFactor, center, stretch)
    image_original = single(image);
    block_size = 2*blockRadius+1;
    mean_val = imfilter(image_original, ones(block_size)/prod(block_size), 'symmetric');
    if stretch
        squared_mean = imfilter(image_original.^2, ones(block_size)/prod(block_size), 'symmetric');
        std_dev = sqrt(squared_mean - mean_val.^2);
        d = meanFactor * std_dev;
    end

    if center && stretch
        min_val = mean_val - d;
        output = ((image_original - min_val) ./ (2 * d) * range(image_original(:)) + min(image_original(:)));
    elseif center
        output = image_original - mean_val + mean(image_original,1:ndims(blockRadius));
    elseif stretch
        output = (image_original - mean_val) / (2 * d) * range(image_original(:)) + mean_val;
    end
end