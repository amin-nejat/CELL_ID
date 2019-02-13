function resid = fit_gaussian_fixed_cov(x, vol, norm, trunc, cov)
% Calculating residual between a multi-color Gaussian function and the
% original image where the covariance of the Gaussian is fixed.
%
% Amin Nejat

    x = x./norm';
    
    mu = x(1:3);
    
    props = x(6:size(vol,4) + 5)*exp(x(4));
    baseline = x(size(vol,4) + 6:end)*x(5);
    sz = [size(vol, 1), size(vol, 2), size(vol, 3), size(vol, 4)];
    
    recon = simulate_gaussian(sz(1: 3), ...
        mu', cov, props', baseline', zeros(1, size(vol, 4)), trunc);
    
    resid = sqrt(sum((vol(:) - recon(:)).^2));
end