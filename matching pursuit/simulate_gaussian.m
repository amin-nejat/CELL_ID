function volume = simulate_gaussian(sz, mu, sigma, props, noise_mean, noise_std, truncate_percent)
% Simulate a truncated Gaussian function for fitting procedure. This
% function is used by matching pursuit algorithm.
%
% Amin Nejat

%% Gaussian function volume
% sz = [25, 25, 25];
% mu = [13, 13, 13];
% sigma = [4, 0, 0; 0, 4, 0; 0, 0, 4];
% props = [0.1, 0.4, 0.1];
% noise_mean = [0.0001, 0.0002, 0.0001];
% noise_std = [0.00002, 0.00001, 0.00001];

[X, Y, Z] = meshgrid(1: sz(1), 1: sz(2), 1: sz(3));

X = permute(X, [2,1,3]);
Y = permute(Y, [2,1,3]);
Z = permute(Z, [2,1,3]);

pos = [X(:) Y(:) Z(:)];

p = mvnpdf(pos, mu, sigma);

volume = zeros(sz);

prob = reshape(p, sz(1: 3));
prob(prob < prctile(prob(:), truncate_percent)) = 0;
prob = prob / max(prob(:));

for ch = 1: length(props)
    volume(:, :, :, ch) = props(ch)*prob + noise_mean(ch) + noise_std(ch)*randn(sz(1: 3));
end



%% Probabilistic volume

% volume = zeros(25, 25, 25);
% X = round(mvnrnd([13, 13, 13], [4, 0, 0; 0, 4, 0; 0, 0, 4], 500000));
% [Xu, ia, ic] = unique(X, 'rows', 'stable');
% h = accumarray(ic, 1);
% 
% ind = sub2ind(size(volume), Xu(:, 1), Xu(:, 2), Xu(:, 3));
% 
% volume(ind) = h;

end