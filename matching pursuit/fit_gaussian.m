function [shape, sp, goodness] = fit_gaussian(bpatch, szext, colorvec, fsize, trunc, absolute_position)
    options = optimoptions('fmincon');
    options.MaxFunEvals = 10000;
    options.MaxIter = 10000;
    options.StepTolerance = 5e-10;
    % options.PlotFcn = 'optimplotx';
    % options.Display = 'iter';

    norm = [10./szext(1: 3), ... % loc
            1, ... % color scale
            1, ... % noise scale
            ones(1, szext(4)), ... % color
            ones(1, szext(4)), ... % noise
            ];

    init_eig = 2.5*sort(fsize)';
    init_bound = init_eig-1;

    nonlcon = @(x) gaussian_confun_cov(x, ones(1,8), init_eig, init_bound);

    colors = reshape(bpatch, [numel(bpatch)/size(bpatch, 4), size(bpatch, 4)]);
    colors(colors == 0) = eps;
    noisevec = prctile(colors, 10);

    colorvec(colorvec < 0) = eps;

    color_level = sum(colorvec(:));
    noise_level = sum(noisevec(:));

    colorvec = colorvec/sum(colorvec(:));
    noisevec = noisevec/sum(noisevec(:));


    x0 =   double([fsize+1, ... % loc
        log(color_level), ... % color scale
        noise_level, ... % noise scale
        colorvec, ... % color
        noisevec, ... % noise
        ]).*norm;


    if isnan(x0(4))
        x0(4) = 2;
    end

    bounds = [30*fsize./(6*szext(1:3)), ... % loc
            2, ... % color scale
            0.2, ... % noise scale
            0.1*ones(1, szext(4)), ... % color
            0.5*ones(1, szext(4)), ... % noise
            ];

    x0(6: 5+2*szext(4)) = max(0, x0(6: 5+2*szext(4)));

    lb = x0 - bounds;
    ub = x0 + bounds;

    lb(6: 5+2*szext(4)) = max(0, lb(6: 5+2*szext(4)));
    lb(4) = max(1, lb(4));

    lb(5) = -1;
    ub(5) = 1;



    A_eq = [zeros(1, 5), ones(1, szext(4)), zeros(1, szext(4)); ...
           [zeros(1, 5), zeros(1, szext(4)), ones(1, szext(4))]];
    b_eq = [1; 1];


    f = @(x) fit_gaussian_fixed_cov(x', bpatch, norm, trunc, eye(3));
    res_fixed_cov = fmincon(f, x0, [], [], A_eq, b_eq, lb, ub, [], options);


    f = @(x) fit_gaussian_cov(x', bpatch, ones(1,8), trunc, res_fixed_cov(1:3)'./norm(1:3)', res_fixed_cov(6:5+szext(4)), res_fixed_cov(6+szext(4):5+2*szext(4)));
    res = fmincon(f, [sqrt([init_eig(3),0,init_eig(2),0,0,init_eig(1)]),res_fixed_cov(4:5)], [], [], [], [], [eps,eps,eps,eps,eps,eps,1,-1], [10,10,10,10,10,10,100,1], nonlcon, options);

    [shape, rec, tr, cov, col] = get_gaussian_cov(res, [2*fsize+1, szext(4)], ones(1,8), trunc, res_fixed_cov(1:3)./norm(1:3), res_fixed_cov(6:5+szext(4)));
    bas = res_fixed_cov(6+szext(4):5+2*szext(4))/res(end);
    
    subplot(1,3,1)
    image(squeeze(max(rec(:,:,:,[4,3,1]), [], 3))/20)
    subplot(1,3,2)
    image(squeeze(max(bpatch(:,:,:,[4,3,1]), [], 3))/20)
    subplot(1,3,3)
    image(squeeze(max(bpatch(:,:,:,[4,3,1]) - rec(:,:,:,[4,3,1]), [], 3))/20)
    drawnow
    
    
    reccoef = corrcoef(rec(:), bpatch(:));
    goodness = reccoef(1,2);
    
    relative_position = (res_fixed_cov(1:3)-x0(1:3))./norm(1:3);
    
    sp.mean = relative_position+absolute_position;
    sp.cov(1,:,:) = cov;
    sp.color = col;
    sp.baseline = bas;
end