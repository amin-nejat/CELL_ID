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

function s = logsumexp(a, dim)
% Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
% Default is dim = 1 (columns).
% logsumexp(a, 2) will sum across rows instead of columns.
% Unlike matlab's "sum", it will not switch the summing direction
% if you provide a row vector.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

	if nargin < 2
	  dim = 1;
	end

	% subtract the largest in each column
	[y, i] = max(a,[],dim);
	dims = ones(1,ndims(a));
	dims(dim) = size(a,dim);
	a = a - repmat(y, dims);
	s = y + log(sum(exp(a),dim));
	i = find(~isfinite(y));
	if ~isempty(i)
	  s(i) = y(i);
	end
end