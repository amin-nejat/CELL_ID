function D=pdist2_maha(X,Y,Sigma)
D=zeros(size(X,1),size(Y,1));
for i=1:size(Y,1)
    D(:,i)=pdist2(X,Y(i,:),'mahalanobis',Sigma(:,:,i));
end