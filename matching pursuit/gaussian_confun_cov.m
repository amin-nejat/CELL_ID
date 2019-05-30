function [c, ceq] = gaussian_confun_cov(x, norm, init, bound)
% Nonlinear inequality constraints (eigenvalues) for fmincon.
%
% Amin Nejat

x = x./norm;

L = zeros(3, 3);
L([1,2,3,5,6,9]) = x(1: 6);
variances = L'*L;

v = sort(eig(variances));

c = -[v'-init+bound, bound-v'+init];

% Nonlinear equality constraints

ceq = [];