function P= sinkhorn_logspace(logP,niter)

a=nanmin(nanmin(logP(~isinf(logP))));
logP(isinf(logP))=a-10000;
n1 = size(logP,1);
n2 = size(logP,2);

logP2 =(a-10000)*ones(max(n1,n2), max(n1,n2));

if(n1<n2)
    logP2(1:n1,:)=logP;
elseif(n1>n2)
    logP2(:,1:n2)=logP;
else
    logP2 = logP;
end

for i=1:niter
        logP2 = logP2 - repmat(logsumexp(logP2, 1), [size(logP2,1), 1]);
        logP2 = logP2 - repmat(logsumexp(logP2, 2), [1, size(logP2,2)]);
end
P = exp(logP2);
P = P(1:n1,1:n2);
end