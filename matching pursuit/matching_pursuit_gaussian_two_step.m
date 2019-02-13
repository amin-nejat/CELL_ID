function supervoxels = matching_pursuit_gaussian_two_step(volume, filter, n_objects, trunc)
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

options = optimoptions('fmincon');
options.MaxFunEvals = 10000;
options.MaxIter = 10000;
options.StepTolerance = 5e-10;
% options.PlotFcn = 'optimplotx';
% options.Display = 'iter';

norm = [10./szext(1: 3), ... % loc
        1, ... % color scale
        1, ... % noise scale
        ones(1, size(volume, 4)), ... % color
        ones(1, size(volume, 4)), ... % noise
        ];
    
init_eig = 2*sort(fsize)';
init_bound = init_eig-1;

nonlcon = @(x) gaussian_confun_cov(x, ones(1,8), init_eig, init_bound);

for object = 1: n_objects
    object
    rho_mag = max(rho, [], 4);
    [~, lmidx] = max(rho_mag(:));
    [x, y, z] = ind2sub(size(rho_mag), lmidx);
    
    bpatch = subcube(volume, [x, y, z], fsize);
    
    colors = reshape(bpatch, [numel(bpatch)/size(bpatch, 4), size(bpatch, 4)]);
    colors(colors == 0) = eps;
    
    colorvec = squeeze(rho(x,y,z,:))';
    noisevec = prctile(colors, 10);
    
    colorvec(colorvec < 0) = eps;
%     noisevec(noisevec < 0) = eps;
    
    color_level = sum(colorvec(:));
    noise_level = sum(noisevec(:));
    
    colorvec = colorvec/sum(colorvec(:));
    noisevec = noisevec/sum(noisevec(:));
    

    x0 =   [fsize+1, ... % loc
        log(color_level), ... % color scale
        noise_level, ... % noise scale
        colorvec, ... % color
        noisevec, ... % noise
        ].*norm;
    
    
    if isnan(x0(4))
        x0(4) = 2;
    end

    bounds = [10*fsize./(6*szext(1:3)), ... % loc
            2, ... % color scale
            0.2, ... % noise scale
            0.1*ones(1, size(volume, 4)), ... % color
            0.5*ones(1, size(volume, 4)), ... % noise
            ];
    
    x0(6: 5+2*size(volume,4)) = max(0, x0(6: 5+2*size(volume,4)));
    
    lb = x0 - bounds;
    ub = x0 + bounds;
    
    lb(6: 5+2*size(volume,4)) = max(0, lb(6: 5+2*size(volume,4)));
    lb(4) = max(1, lb(4));
    
    lb(5) = -1;
    ub(5) = 1;
    
    
    
    A_eq = [zeros(1, 5), ones(1, size(volume, 4)), zeros(1, size(volume, 4)); ...
           [zeros(1, 5), zeros(1, size(volume, 4)), ones(1, size(volume, 4))]];
    b_eq = [1; 1];
    
%     init_eig = 2*sort(fsize)';
%     init_bound = init_eig-1;
    
    
    f = @(x) fit_gaussian_fixed_cov(x', bpatch, norm, trunc, eye(3));
    res_fixed_cov = fmincon(f, x0, [], [], A_eq, b_eq, lb, ub, [], options);
    
    
    f = @(x) fit_gaussian_cov(x', bpatch, ones(1,8), trunc, res_fixed_cov(1:3)'./norm(1:3)', res_fixed_cov(6:5+size(volume,4)), res_fixed_cov(6+size(volume,4):5+2*size(volume,4)));
    res = fmincon(f, [sqrt([init_eig(3),0,init_eig(2),0,0,init_eig(1)]),res_fixed_cov(4:5)], [], [], [], [], [eps,eps,eps,eps,eps,eps,1,-1], [10,10,10,10,10,10,100,1], nonlcon, options);
    
    [shape, rec, tr, cov, col] = get_gaussian_cov(res, [2*fsize+1, szext(4)], ones(1,8), trunc, res_fixed_cov(1:3)./norm(1:3), res_fixed_cov(6:5+size(volume,4)));
    bas = res_fixed_cov(6+size(volume,4):5+2*size(volume,4))/res(end);
    
%     imshow3D(cat(2, bpatch(:,:,:,[1,2,3]), ...
%         rec(:,:,:,[1,2,3]), ...
%         bpatch(:,:,:,[1,2,3]) - rec(:,:,:,[1,2,3]))/20);
%     drawnow;
    
    fshape = imfilter(shape, filter, 'full');
    residual = placement(szext(1:3), [x,y,z], fshape);
    
    for ch = 1: size(rho, 4)
        rho(:, :, :, ch) = rho(:, :, :, ch)-residual*col(ch);
    end
    
    supervoxels.mean(object, :) = (res_fixed_cov(1:3)-x0(1:3))./norm(1:3)+[x,y,z];
    supervoxels.cov(object, :, :) = cov;
    supervoxels.color(object, :) = col;
    supervoxels.baseline(object, :) = bas;
    
end



end

