function [shape, recon, mu, cov, colors] = get_gaussian_cov(x, sz, norm, trunc, mu, color_props)
% Creating Gaussian function for a given set of parameters, used for
% visualization.
%
% Amin Nejat

    x = x./norm;
    
    L = zeros(3, 3);
    L([1,2,3,5,6,9]) = x(1:6);
    cov = L*L';
    
    colors = color_props*exp(x(7));
    
    recon = simulate_gaussian(sz, ...
        mu, cov, colors, zeros(1, sz(4)), zeros(1, sz(4)), trunc);

    shape = recon(:, :, :, 1)/colors(1);
end