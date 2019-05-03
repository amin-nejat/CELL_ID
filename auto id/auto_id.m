function auto_id(sp, method)
% Automatic labeling of previously identified neurons based on mixture of
% Gaussian model.
%
% sp is supervoxel data structure that is output by MP algorithm
%
% Erdem

compartment = lower(sp.bodypart);
positions=sp.get_positions();
colors = sp.get_colors();

if strcmp(lower(method), 'mog')
    
    
    %% MoG Based
    D=squareform(pdist(bsxfun(@times,positions,[1 1 4])));
    
    for j=1:size(colors, 1)
        worm_colors{j}=colors(:,[1 2 3]);
        worm_distances{j}=D(j,:);
        worm_positions{j}=positions(j,:);
    end
    
    
    %% Training ID likelihoods
    load(['mog_model_', compartment, '.mat']);
    neurons=N;
    N=length(worm_distances);
    for k=1:N
        k
        for j=1:length(mog)
            try
                LL(k,j)=nansum(log(pdf(mog{j},[worm_distances{k}((1:N)~=k)' worm_colors{k}((1:N)~=k,:)]))) + log(mvnpdf(worm_colors{k}(k,:),cd{j}.mean,cd{j}.cov));
            catch
            end
        end
    end
    
    
    
    %% Linear assignment
    a=nanmin(nanmin(LL(~isinf(LL))));
    LL(isinf(LL))=a-10000;
    
elseif strcmp(lower(method),'rwc_old')
    
    %% RWC Based
    
    load(['aligned_worms_', compartment, '.mat']); %% LOAD ALIGNED TRAINING DATA HEAD
    neurons = N;
    sigma=inv(nancov(reshape(permute(M,[1 3 2]),[size(M,1)*size(M,3) size(M,2)])));
    
    X=[positions.*[1 1 4] colors(:,1:3)]; % Make sure that z-pixels are scaled accordingly
    Phat=zeros(size(M,1),size(X,1));
    for j=1:size(M,3)
        Y=M(:,:,j);
        neurons=N;
        idx=find(~any(isnan(Y),2));
        neurons(any(isnan(Y),2))=[];
        Y(any(isnan(Y),2),:)=[];
        [beta,P,Xhat,cost]=regression_woc(X,Y,sigma,20);
        Phat(idx,:)=Phat(idx,:)+P;
        Xmatch=P*Xhat;
        
        
        
        
    end
    
    LL=-pdist2(Xhat,nanmean(M,3),'mahalanobis',inv(sigma)).^2;
    
    a=nanmin(nanmin(LL(~isinf(LL))));
    LL(isinf(LL))=a-10000;
    
elseif strcmp(lower(method),'rwc')
    
    load(['model3_', compartment, '.mat']); %% LOAD ALIGNED TRAINING DATA HEAD
    colors=colors(:,[1 2 3]);
    positions=positions(:,[1 2 3]);
    neurons = N;
    iter=3;
    %% Anonymous functions
    vec = @(x)(x(:));
    sum_square = @(x)(sum(vec(x.^2)));
    
    
    %% initialization
    tmp.M=model.M0;
    [~,LL2]=model_2_test(colors,positions,tmp);
    P=munkres(-LL2'-logsumexp(LL2',2));
    M=model.M;
    Sigma=model.Sigma;
    
    X=[positions colors];
    X(any(isnan(X),2),:)=[];
    X0=X;
    Y=nanmean(M,3);
    beta = linsolve(P'*[X ones(size(X,1),1)],Y);
    X = [X ones(size(X,1),1)]*beta;
    LLhat=zeros(size(M,1),size(X,1));
    for j=1:size(M,3)
        Y=M(:,:,j);
        idx=find(~any(isnan(Y),2));
        Y(any(isnan(Y),2),:)=[];
        beta=[eye(size(X,2));zeros(1,size(Y,2))];
        for t=1:iter
            X = [X ones(size(X,1),1)]*beta;
            positions = X(:,1:size(positions,2)); colors=X(:,size(positions,2)+1:end);
            [P2,~]=model_2_test(colors,positions,model);
            P=P2(idx,:)';
            beta = linsolve(P'*[X ones(size(X,1),1)],Y);
        end
        Xhat=[X ones(size(X,1),1)]*beta;
        LL=-pdist2_maha(Xhat,Y,Sigma).^2;
        LLhat(idx,:)=LLhat(idx,:)+LL';
    end
    LL=LLhat-logsumexp(LLhat,1);
end


sp.add_meta_data('id_method', method);
sp.add_meta_data('LL', LL);
sp.add_meta_data('ids', neurons);

end