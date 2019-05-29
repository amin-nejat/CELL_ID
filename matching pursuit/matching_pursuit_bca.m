function sp = matching_pursuit_bca(volume, sp, fsize, trunc, iterations)
% Block coordinate ascent updates of the MoG parameters found by MP.
%
% Amin Nejat

resid = double(volume);
n_objects = size(sp.mean,1);
recon = cell(n_objects,1);

szext = [size(volume), 1];

options = optimoptions('fmincon');
options.MaxFunEvals = 1000;
options.MaxIter = 1000;
options.StepTolerance = 5e-10;

norm = [10./szext(1: 3), ... % loc
        1, ... % color scale
        1, ... % noise scale
        ones(1, size(volume, 4)), ... % color
        ones(1, size(volume, 4)), ... % noise
        ];
    
init_eig = sort(fsize);
init_bound = init_eig-1;

nonlcon = @(x) gaussian_confun_cov(x, ones(1,8), init_eig, init_bound);


for object=1:n_objects
    rec = simulate_gaussian(2*fsize+1, ... % size
            fsize+1+sp.mean(object,:)-round(sp.mean(object,:)), ... % center
            squeeze(sp.cov(object,:,:)), ... % cov
            sp.color(object,:), ... % colors
            0*sp.baseline(object,:), ... % baseline mean
            0*sp.baseline(object,:), ... % baseline std
            trunc); % truncation
        
    recon{object} = rec;
    resid = resid - placement(szext(1:3), round(sp.mean(object,:)), recon{object});
end


for iter = 1: iterations
    permutation = randperm(n_objects);
    for index = 1: length(permutation)
        object = permutation(index);
        mean = sp.mean(object, :);
        resid = resid + placement(szext(1:3), round(mean), recon{object});
        bpatch = subcube(resid, round(mean), fsize);

    %     imshow3D(resid(:,:,:,[1,2,2],:))
    %     drawnow;

        colorvec = sp.color(object, :);
        noisevec = sp.baseline(object, :);

        color_level = sum(colorvec);
        noise_level = sum(noisevec);

        colorvec = colorvec/sum(colorvec(:));
        noisevec = noisevec/sum(noisevec(:));


        x0 =   [mean-round(mean)+fsize+1, ... % loc
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

        f = @(x) fit_gaussian_fixed_cov(x', bpatch, norm, trunc, squeeze(sp.cov(object,:,:)));
        res_fixed_cov = fmincon(f, x0, [], [], A_eq, b_eq, lb, ub, [], options);

        cholesky = chol(squeeze(sp.cov(object,:,:)));

        f = @(x) fit_gaussian_cov(x', bpatch, ones(1,8), trunc, res_fixed_cov(1:3)'./norm(1:3)', res_fixed_cov(6:5+size(volume,4)), res_fixed_cov(6+size(volume,4):5+2*size(volume,4)));
        res = fmincon(f, [cholesky([1, 4, 5, 7, 8, 9]),res_fixed_cov(4:5)], [], [], [], [], [0,0,0,0,0,0,1,-1], [10,10,10,10,10,10,100,1], nonlcon, options);

        [shape, rec, tr, cov, col] = get_gaussian_cov(res, [2*fsize+1, szext(4)], ones(1,8), trunc, res_fixed_cov(1:3)./norm(1:3), res_fixed_cov(6:5+size(volume,4)));
        bas = res_fixed_cov(6+size(volume,4):5+2*size(volume,4))/res(end);



    %     imshow3D(cat(2, bpatch(:,:,:,[1,1,2]), ...
    %         rec(:,:,:,[1,1,2]), ...
    %         bpatch(:,:,:,[1,1,2])-rec(:,:,:,[1,1,2]))/20);
    %     drawnow;




        sp.mean(object, :) = (res_fixed_cov(1:3)-x0(1:3))./norm(1:3)+round(mean);
        sp.cov(object, :, :) = cov;
        sp.color(object, :) = col;
        sp.baseline(object, :) = bas;

        recon{object} = rec;
        resid = resid - placement(szext(1:3), sp.mean(object, :), rec);
    end
end

end