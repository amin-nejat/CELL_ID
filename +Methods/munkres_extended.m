function [P, cost, P_full] = munkres_extended(Cost, p, n1)

% Solves an assignment problem, but accounts for possible n1 false
% positives, so each detect is a false positive with probability p/n1.
n2 = size(Cost,2);
p2 = p/n1;
p1 = (1-p)/n2;
Cost = Cost - log(p1);


Cost(:,n2+1:n2+n1) = -log(p2);


[P, cost] = Methods.munkres(Cost);
P_full = P;
P= P_full(:,1:n2);

