function [shape, recon, mu, cov, colors, baseline] = get_gaussian_fixed_cov(x, sz, rel_mu, norm, trunc, cov)
% Creating Gaussian function for a given set of parameters where the 
% covariance is fixed, used for visualization.
%
% Amin Nejat

    x = x./norm;
    
    mu = x(1: 3)+rel_mu;
    
    colors = x(6: sz(4) + 5)*exp(x(4));
    baseline = x(sz(4) + 6: end)*x(5);
    
    recon = simulate_gaussian(sz, ...
        mu, cov, colors, zeros(1, sz(4)), zeros(1, sz(4)), trunc);

    shape = recon(:, :, :, 1)/colors(1);
end