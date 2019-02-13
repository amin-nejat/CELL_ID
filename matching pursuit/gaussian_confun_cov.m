function [c, ceq] = gaussian_confun_cov(x, norm, init, bound)
% Nonlinear inequality constraints (eigenvalues) for fmincon.
%
% Amin Nejat

x = x./norm;

L = zeros(3, 3);
L([1, 4, 5, 7, 8, 9]) = x(1: 6);
L = reshape(L, 3, 3);
variances = L'*L;

v = sort(eig(variances));

c = -[v'-init+bound, bound-v'+init];

% Nonlinear equality constraints

ceq = [];