function supervoxels = matching_pursuit_gaussian_two_step(volume, filter, n_objects, trunc, cov_threshold)
% Matching pursuit algorithm for finding the best mixture of Gaussians that
% fits to input dataset. The algorithm runs greedy by subtracting off the
% brightest Gaussian at each iteration.
%
% Amin Nejat

supervoxels = [];

volume = double(volume);

szext = [size(volume), 1];
fsize = size(filter)-1;

rho = filter_frame(volume, filter);

N = 0;

while N < n_objects && max(rho(:)) > 0.1
    disp(['Object: ', num2str(N)]);
    max(rho(:)) 
    
    rho_mag = max(rho, [], 4);
    [~, lmidx] = max(rho_mag(:));
    [x,y,z] = ind2sub(size(rho_mag), lmidx);
    
    bpatch = subcube(volume, [x,y,z], fsize);
    
    [shape, sp, goodness] = fit_gaussian(bpatch, szext, squeeze(rho(x,y,z,:))', fsize, trunc, [x,y,z]);
    
    fshape = imfilter(shape, filter, 'full');
    residual = placement(szext(1:3), [x,y,z], fshape);
    
    for ch = 1: size(rho, 4)
        rho(:, :, :, ch) = rho(:, :, :, ch)-residual*sp.color(ch);
    end
     
    
    if max(eig(squeeze(sp.cov))) > cov_threshold
        supervoxels = union_sp(supervoxels, sp);
        N = size(supervoxels.mean, 1);
    end
end



end

