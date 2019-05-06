function beta=MCR_solver(Y,X,sigma)
% MCR - Multiple covariance regression solver
% Y - Target (n x d)
% X - Source (n x p)
% Sigma - covariances for each row of Y (d x d x n)
% Solving: \sum_i \|Y_i - X_i \beta\|_{sigma_i}^2

for i=1:size(Y,1)
    A(:,:,i)=kron(inv(sigma(:,:,i)),X(i,:)'*X(i,:));
    B(:,:,i)=X(i,:)'*Y(i,:)*inv(sigma(:,:,i));
end

beta=reshape(inv(nansum(A,3))*vec(nansum(B,3)),[size(X,2),size(Y,2)]);
end