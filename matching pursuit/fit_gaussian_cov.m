function resid = fit_gaussian_cov(x, vol, norm, trunc, mu, color_props, baseline_props)
% Calculating residual between a multi-color Gaussian for a given set of
% parameters and the original image. In this function covariance of the
% Gaussian is the only variable, rest of the parameters are fixed. Used by 
% fmincon to find the best fit Gaussian function.
%
% Amin Nejat
    
    x = x./norm';
    
    L = zeros(3, 3);
    L([1, 4, 5, 7, 8, 9]) = x(1:6);
    variances = L*L';
    
    props = color_props*exp(x(7));
    baseline = baseline_props*x(8);
    sz = [size(vol, 1), size(vol, 2), size(vol, 3), size(vol, 4)];
    
    recon = simulate_gaussian(sz(1: 3), ...
        mu', variances, props', baseline', zeros(1, size(vol, 4)), trunc);
    
    resid = sqrt(sum((vol(:) - recon(:)).^2));
end