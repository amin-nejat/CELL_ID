function [beta,P,Xhat,cost]=regression_woc(X,Y,sigma,iter)

beta=[eye(size(X,2));zeros(1,size(Y,2))];
for t=1:iter
    X = [X ones(size(X,1),1)]*beta;
    D=pdist2(Y,X,'mahalanobis',inv(sigma));
    [P,cost(t)]=munkres(D);
    beta = linsolve(P*[X ones(size(X,1),1)],Y);
end
Xhat=X;
end