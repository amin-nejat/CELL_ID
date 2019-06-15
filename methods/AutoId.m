classdef AutoId < handle
    %AUTO_ID Auto.
    %
    %   An auto_id object is has a set of possible names (different if tail/head), the method used, a likelihood
    %   matrix, and a set of ids for each of observed neurons. 

    

    properties
        neuron_names %set of possible neuron names, depending on head or tail. 
        id_method %mog or rwoc
        LL %likelihood matrix, of size n_mp x n_possible
        ids %list of strings of inferred neuron names of size n_mp
        ids_prob %array of strings of most likely neuron names of size n_mp x 7 (7 most likely)
        assignments %matrix of assignment for each observed neuron to canonical identity. Binary matrix of size n_mp x n_possible
        assignment_det %vector of size n_obs such that assignment_det(i)=j for the entry such that assignments(i,j)=1;
        assignment_prob_ranks %n_obsx7 matrix with the 7 most likely assignment
        assignment_prob_probs %n_obsx7 matrix with the corresponding probabilities of assignment_prob_ranks
    end
    
    methods(Static)
        
        function obj = instance(image, type)
             persistent instance
             %if isempty(instance)
                obj = AutoId(image, type);
              %  instance = obj;
             %else
               % obj = instance;
             %end
        end
      
        function D=pdist2_maha(X,Y,Sigma)
            D=zeros(size(X,1),size(Y,1));
            for i=1:size(Y,1)
                D(:,i)=pdist2(X,Y(i,:),'mahalanobis',Sigma(:,:,i));
            end
        end
        
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
        
        function [S,R,T]=scaled_rotation(X,Y)
        %Solves for Y = S*R*X + T
        %where S is a diagonal scaling matrix
        %R is a rotation matrix i.e. orthonormal and det(R)=1
        %T is a translation

            % Remove NAN rows
            idx=find(all(~isnan([X Y]),2));
            X=X(idx,:);
            Y=Y(idx,:);


            % De-mean
            Yhat=Y-mean(Y,1);
            Xhat=X-mean(X,1);

            % Scale
            sx = sqrt(sum(sum(Xhat.^2,2),1)/size(Xhat,1));
            sy = sqrt(sum(sum(Yhat.^2,2),1)/size(Yhat,1));

            Yhat=Yhat./sy;
            Xhat=Xhat./sx;

            % Solve rotation
            C = Yhat'*Xhat;
            [U,~,V]=svd(C);

            R0=V*U';


            % Put it all together

            S=diag(sy./sx);
            R=R0;
            T=(mean(Y,1)-mean(X,1)*R0*S);
        end
        
        function P= sinkhorn_logspace(logP,niter)
            a=nanmin(nanmin(logP(~isinf(logP))));
            logP(isinf(logP))=a-10000;
            logP(isnan(logP))=-a-10000;
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
    
    
    
        function [neuron_names, LL, id_method] = auto_id(sp, method)
        % Automatic labeling of previously identified neurons based on mixture of
        % Gaussian model.
        %
        % sp is supervoxel data structure that is output by MP algorithm
        %
        % Erdem

            compartment = lower(sp.bodypart);
            positions=sp.get_positions();
            colors = sp.get_colors();
            [labels, conf] = sp.get_human_labels();
            
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
                    [beta,P,Xhat,cost] = AutoId.regression_woc(X,Y,sigma,20);
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
                known = AutoId.find_labeled_index(labels,conf, neurons);
                iter=2;
                %% Anonymous functions
                vec = @(x)(x(:));
                sum_square = @(x)(sum(vec(x.^2)));
                
                
                %% initialization
                
                n_total_iter =100;
                M=model.M;
                M0=model.M0;
                sigma0 = nancov(reshape(permute(M,[1 3 2]),[size(M,1)*size(M,3) size(M,2)]));
                sigma00 = nancov(reshape(permute(M0,[1 3 2]),[size(M0,1)*size(M0,3) size(M0,2)]));
                
                
                Sigma_eval = model.Sigma;
                lambda=0.5;
                for i=1:size(M,1)
                    
                    Sigma(:,:,i)=lambda*model.Sigma(:,:,i)+(1-lambda)*sigma0;
                    Sigma0(:,:,i)=sigma0;
                    Sigma00(:,:,i)=sigma00;
                    
                end
                
                model.Sigma = Sigma;
                tmp.M = model.M0;
                tmp.Sigma  = Sigma00;
                
                
                P=AutoId.model_2_test(colors,positions,tmp,known);
                
                X=[positions.*[1 1 4] colors(:,1:3)]; % Make sure that z-pixels are scaled accordingly
                
                X(any(isnan(X),2),:)=[];
                
                YM =nanmean(M,3);
                beta = linsolve(P'*[X ones(size(X,1),1)],YM);
                X = [X ones(size(X,1),1)]*beta;
                
                
                LLhat=zeros(size(M,1),size(X,1));
                
                
                for jj=1:n_total_iter
                    
                    disp(['Iteration ' num2str(jj) '/' num2str(n_total_iter) ' done']);
                %                     
                    j = randsample(size(M,3),1);
                    Y = M(:,:,j);
                    idx = find(~any(isnan(Y),2));
                    Y(any(isnan(Y),2),:)=[];
                    
                    
                    LL =  -pdist2_maha(X,YM,Sigma_eval);
                    LLhat(idx,:)=LL';
                    
                    
                    
                    beta=[eye(size(X,2));zeros(1,size(Y,2))];
                    
                    for t=1:iter
                        X = [X ones(size(X,1),1)]*beta;
                        positions = X(:,1:size(positions,2));
                        colors = X(:,size(positions,2)+1:end);
                        P = AutoId.model_2_test(colors,positions,model,known);
                        beta = linsolve(P(:,idx)'*[X ones(size(X,1),1)],Y);
                        
                    end
                    
                    
                end
                
                %                 tmp.M=model.M0;
                %                 [~,LL2]=AutoId.model_2_test(colors,positions,tmp);
                %                 P=munkres(-LL2'-logsumexp(LL2',2));
                %                 M=model.M;
                %                 Sigma=model.Sigma;
                %
                %                 X=[positions colors];
                %                 X(any(isnan(X),2),:)=[];
                %                 X0=X;
                %                 Y=nanmean(M,3);
                %                 beta = linsolve(P'*[X ones(size(X,1),1)],Y);
                %                 X = [X ones(size(X,1),1)]*beta;
                %                 LLhat=zeros(size(M,1),size(X,1));
                %                 for j=1:size(M,3)
                %                     disp([num2str(j) '/' num2str(size(M,3)) ' done']);
                %                     Y=M(:,:,j);
                %                     idx=find(~any(isnan(Y),2));
                %                     Y(any(isnan(Y),2),:)=[];
                %                     beta=[eye(size(X,2));zeros(1,size(Y,2))];
                %                     for t=1:iter
                %                         X = [X ones(size(X,1),1)]*beta;
                %                         positions = X(:,1:size(positions,2)); colors=X(:,size(positions,2)+1:end);
                %                         [P2,~]=AutoId.model_2_test(colors,positions,model);
                %                         P=P2(idx,:)';
                %                         beta = linsolve(P'*[X ones(size(X,1),1)],Y);
                %                     end
                %                     Xhat=[X ones(size(X,1),1)]*beta;
                %                     LL=-AutoId.pdist2_maha(Xhat,Y,Sigma).^2;
                %                     LLhat(idx,:)=LLhat(idx,:)+LL';
                %                 end
                %                 LL=(LLhat-logsumexp(LLhat,1))';
                
                
                
            end
            
            
            id_method = method;
            LL =  LL;
            neuron_names = neurons;
            
        end
        function [P] =model_2_test(colors,positions,model, known)
            
            M=model.M;
            X=[positions colors];
            X(any(isnan(X),2),:)=[];
            
            LLhat=zeros(size(M,1),size(X,1));
            order = randperm(size(M,3));
            for jj=1:length(order)
                j=order(jj);
                Y=M(:,:,j);
                idx=find(~any(isnan(Y),2));
                Y(any(isnan(Y),2),:)=[];
                D = pdist2_maha(X,Y,model.Sigma);
                LLhat(idx,:)=LLhat(idx,:)+D';
                
            end
            LLhat = LLhat';
            if(~isempty(known))
                LLhat(known(:,1),:) = 10000000000;
                for i=1:length(known)
                    LLhat(known(i,1),known(i,2))=-100000000;
                end
            end
            P = munkres(LLhat);
        end
%         function [assignments,LL,beta]=model_2_test(colors,positions,model)
%             M=model.M;
%             X=[positions colors];
%             X(any(isnan(X),2),:)=[];
%             sigma0=inv(nancov(reshape(permute(M,[1 3 2]),[size(M,1)*size(M,3) size(M,2)])));
%             Phat=zeros(size(M,1),size(X,1));
%             for j=1:size(M,3)
%                 Y=M(:,:,j);
%                 idx=find(~any(isnan(Y),2));
%                 Y(any(isnan(Y),2),:)=[];
%                 [beta,P,~,~]=AutoId.regression_woc(X,Y,sigma0,3);
%                 Phat(idx,:)=Phat(idx,:)+P;
%             end
% 
%             LL=log(Phat);
%             [assignments]=munkres(-(LL-logsumexp(LL,1)));
%         end

 function known = find_labeled_index(labels, conf, neurons)
            cont=1;
            
            if(~isempty(conf))
            for i=1:length(conf)
                if(conf{i}>0.5)
                name = labels{i};
                known(cont,1)=i;
                known(cont,2)=find(strcmp(name,neurons));
                cont=cont+1;
                end
            end
            else
                known=[];
            end
            if(cont==1)
                known=[];
            end
        end
    end
        
    
    methods
        function obj = AutoId(image, type)
            %class constructor, defines main properties of Auto_Id
            [neuron_names, LL, id_method] = AutoId.auto_id(image, type); %call the main engine of auto_id
            %define main properties of auto_id
            obj.neuron_names = neuron_names; 
            obj.LL = LL;
            obj.id_method = id_method;
           
        end
        
        function value = get_meta_data(obj, key)
            %GET_META_DATA returns the value paired with key.
            value = obj.meta_data(key);
        end
        
      
        function image = add_to_image(obj,image)
            % add_to_image  function simply uptdates some properties of the image
            % structure that depend on auto_id

            image.add_deterministic_ids(obj.ids); %
            image.add_probabilistic_ids(obj.ids_prob);
            image.add_probabilistic_probs(obj.assignment_prob_probs)

            [~,rank_uncertain]=sort(obj.assignment_prob_probs(:,1));
            rank_uncertain_inv(rank_uncertain) = 1:length(rank_uncertain);
            image.add_ranks(rank_uncertain_inv');

        end
        
        function find_ids(obj)
            % find_ids computes some secondary properties of auto_id object
            assignment_prob_ranks = obj.assignment_prob_ranks;
            assignments = obj.assignments;
            neurons = obj.neuron_names;

            ncandidates = size(assignment_prob_ranks, 2);
            ids = repmat({'NaN'}, size(assignments,1),1);
            ids_prob = repmat({'NaN'}, size(assignments,1),ncandidates);

            [rows,cols] = find(assignments);
            ids(rows) = neurons(cols);

            for i=1:size(assignment_prob_ranks,1)
                for j=1:ncandidates
                    if(assignment_prob_ranks(i,j)>0)

                        ids_prob{i,j}=neurons{assignment_prob_ranks(i,j)};
                    end
                end
            end

            ids(cellfun(@isempty, ids)) = {'NaN'};
            ids_prob(cellfun(@isempty, ids_prob)) = {'NaN'};

            obj.ids = ids;
            obj.ids_prob = ids_prob;

        end
       
        function compute_assignments(obj)
           % compute_assignments computes some secondary properties of auto_id object
            LL = obj.LL;
            
            [assignments,~]=munkres(-LL);
            
            ncandidates = 7;
            niter=2000;
            P = AutoId.sinkhorn_logspace(LL, niter);
            
            
            for i=1:size(P,1);
                prow = [P(i,:) 1-nansum(P(i,:))];
                total_prob(i)=nansum(P(i,:));
                [sort_prob,ind_prob] = sort(-prow);
                ind_prob(ind_prob==length(prow))=-1;
                
                assignment_prob_ranks(i,:)=[ind_prob(1:ncandidates)];
                assignment_prob_probs(i,:)=[-sort_prob(1:ncandidates)];
                
                det = find(assignments(i,:));
                
                if(isempty(det))
                    assignment_det(i)=-1;
                else
                    assignment_det(i)=det;
                    
                end
            end
            for i=1:size(P,1)
                ind_det_in_probs = find(assignment_det(i)==assignment_prob_ranks(i,:));
                rest_indices = setdiff([1:ncandidates], ind_det_in_probs);
                rest_indices = rest_indices(1:ncandidates-1);
                ranks = assignment_prob_ranks(i,:);
                probs = assignment_prob_probs(i,:);
                
                assignment_prob_ranks(i, 1) = assignment_det(i);
                assignment_prob_ranks(i, 2:ncandidates)= ranks(rest_indices);
                
            end

            [~,rank_uncertain]=sort(assignment_prob_probs(:,1));
            obj.assignments = assignments;
            obj.assignment_det = assignment_det';
            obj.assignment_prob_ranks = assignment_prob_ranks;
            obj.assignment_prob_probs = assignment_prob_probs;
            
        end
        
        
        function update_auto_id_by_human(obj,image, neuron_i, neuron)
                % update_auto_id_by_human changes the properties of auto_id given human has given 
                % the identity neuron to the neuron_i-th neuron (with the mp order).
                %it also updates the image object accordingly.
                LL = obj.LL;
                names = obj.neuron_names;
                LL(neuron_i,:)= -1000000000;
                
                if ~strcmp(neuron, 'NaN')
                    LL(neuron_i, strcmp(names, neuron))  = 100000;
                end
                
                obj.LL = LL;
                
                obj.compute_assignments();
                obj.find_ids();
                image = add_to_image(obj, image);
        end
        
    end
end

