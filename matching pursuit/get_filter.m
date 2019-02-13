function filter = get_filter(sz, sigma_factor, truncate_percentage)
% Creating a Gaussian filter 2D, 3D, or 4D for smoothing an image.
%
% Amin Nejat

    fmu = (sz+1)/2;
    fsigma = diag(sz)./sigma_factor;
    
    pos = [];
    
    if length(sz) == 2
        [pos(:, 1), pos(:, 2)] = ind2sub(sz, find(ones(sz)));
    elseif length(sz) == 3
        [pos(:, 1), pos(:, 2), pos(:, 3)] = ind2sub(sz, find(ones(sz)));
    elseif length(sz) == 4
        [pos(:, 1), pos(:, 2), pos(:, 3), pos(:, 4)] = ind2sub(sz, find(ones(sz)));
    end
    
    filt = reshape(mvnpdf(pos, fmu, fsigma), sz);
    filt(filt < prctile(filt(:), truncate_percentage)) = 0;

    filter = filt/sqrt(sum(filt(:).*filt(:)));
end