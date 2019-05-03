function [assignments,LL,beta]=model_2_test(colors,positions,model)

M=model.M;
X=[positions colors];
X(any(isnan(X),2),:)=[];
sigma0=inv(nancov(reshape(permute(M,[1 3 2]),[size(M,1)*size(M,3) size(M,2)])));
% for i=1:size(M,1)
%     sigma(:,:,i)=blkdiag(inv(nancov(squeeze(M(i,[1 2 3],:))')),inv(nancov(squeeze(M(i,[4 5 6],:))')));
%     [a,b]=chol(sigma(:,:,i));
%     if b>0
%         sigma(:,:,i)=sigma0;
%     end
% end
Phat=zeros(size(M,1),size(X,1));
for j=1:size(M,3)
    Y=M(:,:,j);
    idx=find(~any(isnan(Y),2));
    Y(any(isnan(Y),2),:)=[];
    [beta,P,~,~]=regression_woc(X,Y,sigma0,3);
    Phat(idx,:)=Phat(idx,:)+P;
end

LL=log(Phat);
% [assignments]=munkres(-(Phat./sum(Phat,1)));
[assignments]=munkres(-(LL-logsumexp(LL,1)));
end