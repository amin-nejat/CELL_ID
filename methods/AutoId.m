classdef AutoId < handle
    %AUTO_ID Auto.
    %
    %   An auto_id object is has a set of possible names (different if tail/head), the method used, a likelihood
    %   matrix, and a set of ids for each of observed neurons. 
    properties(Constant)
        % constants
        ncandidates = 7;
        annotation_weight = 10;
        
        iter_sinkhorn   = 2000;
        iter_regression = 2;
        iter_rwc        = 100;
        
        lambda          =0.5;
        
        min_log_likelihood          = -1e9;
        max_log_likelihood          = 1e5;

    end
    
    properties
        
        atlas = getfield(load('atlas.mat'), 'atlas')
        
        log_likelihood % likelihood matrix, of size n_mp x n_possible
        assignments % matrix of assignment for each observed neuron to canonical identity. Binary matrix of size n_mp x n_possible
        assignment_prob_ranks % n_obsx7 matrix with the 7 most likely assignment
        assignment_prob_probs % n_obsx7 matrix with the corresponding probabilities of assignment_prob_ranks
    end
    
    methods(Static)
        
        function obj = instance()
             persistent instance
             if isempty(instance)
                obj = AutoId();
               instance = obj;
             else
               obj = instance;
             end
        end
      
        function D = pdist2_maha(X,Y,Sigma)
            D=zeros(size(X,1),size(Y,1));
            for i=1:size(Y,1)
                D(:,i)=pdist2(X,Y(i,:),'mahalanobis',Sigma(:,:,i));
            end
        end
        
        function [beta,P,Xhat]=regression_woc(X,Y,sigma,iter)
            beta=[eye(size(X,2));zeros(1,size(Y,2))];
            for t=1:iter
                X = [X ones(size(X,1),1)]*beta;
                D=pdist2(Y,X,'mahalanobis',inv(sigma));
                [P,~]=munkres(D);
                beta = linsolve(P*[X ones(size(X,1),1)],Y);
            end
            Xhat=X;
        end
        
        function vec_x = vec(x)
            vec_x = x(:);
        end
        
        function beta = MCR_solver(Y,X,sigma)
        % MCR - Multiple covariance regression solver
        % Y - Target (n x d)
        % X - Source (n x p)
        % Sigma - covariances for each row of Y (d x d x n)
        % Solving: \sum_i \|Y_i - X_i \beta\|_{sigma_i}^2
            d1 = size(Y,2);
            d2 = size(X,2);
            n = size(Y,1);
            
            A = zeros(d1*d2,d1*d2,n);
            B = zeros(d2,d1,n);
            
            for i=1:size(Y,1)
                A(:,:,i)=kron(inv(sigma(:,:,i)),X(i,:)'*X(i,:));
                B(:,:,i)=X(i,:)'*Y(i,:)/sigma(:,:,i);
            end
            beta=reshape(nansum(A,3)\AutoId.vec(nansum(B,3)),[size(X,2),size(Y,2)]);
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
        
        function P = sinkhorn_logspace(logP,niter)
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
        
        function P = update_permutation(colors, positions, model, known)
            M = model.M;
            X=[positions colors];
            X(any(isnan(X),2),:)=[];
            
            log_likelihoodhat=zeros(size(M,1),size(X,1));
            order = randperm(size(M,3));
            for jj=1:length(order)
                j=order(jj);
                Y=M(:,:,j);
                idx=find(~any(isnan(Y),2));
                Y(any(isnan(Y),2),:)=[];
                D = AutoId.pdist2_maha(X,Y,model.Sigma);
                log_likelihoodhat(idx,:)=log_likelihoodhat(idx,:)+D';
                
            end
            log_likelihoodhat = log_likelihoodhat';
            if(~isempty(known))
                log_likelihoodhat(known(:,1),:) = AutoId.max_log_likelihood;
                for i=1:size(known,1)
                    log_likelihoodhat(known(i,1),known(i,2)) = AutoId.min_log_likelihood;
                end
            end
            P = munkres(log_likelihoodhat);
        end
        
    end
        
    
    methods
        function obj = AutoId()
            % load the statistical atlas
        end
        
        function value = get_meta_data(obj, key)
            %GET_META_DATA returns the value paired with key.
            value = obj.meta_data(key);
        end
        
      
        function add_to_image(obj, im)
            % add_to_image  function simply uptdates some properties of the image
            % structure that depend on auto_id
            
            neuron_names = obj.atlas.(lower(im.bodypart)).N;
            % convert assignment (probabilistic/deterministic) to neuron
            % names
            ids = repmat({'NaN'}, size(obj.assignments,1),1);
            ids_prob = repmat({'NaN'}, size(obj.assignments,1),obj.ncandidates);

            [rows, cols] = find(obj.assignments);
            ids(rows) = neuron_names(cols);

            for i=1:size(obj.assignment_prob_ranks,1)
                for j=1:obj.ncandidates
                    if(obj.assignment_prob_ranks(i,j)>0)
                        ids_prob{i,j} = neuron_names{obj.assignment_prob_ranks(i,j)};
                    end
                end
            end

            ids(cellfun(@isempty, ids)) = {'NaN'};
            ids_prob(cellfun(@isempty, ids_prob)) = {'NaN'};
            
            % update image information
            im.add_deterministic_ids(ids);
            im.add_probabilistic_ids(ids_prob);
            im.add_probabilistic_probs(obj.assignment_prob_probs)

            [~,rank_uncertain] = sort(obj.assignment_prob_probs(:,1));
            rank_uncertain_inv(rank_uncertain) = 1:length(rank_uncertain);
            im.add_ranks(rank_uncertain_inv');
        end
        
       
        function compute_assignments(obj)
           % compute_assignments computes some secondary properties of auto_id object
            [obj.assignments,~] = munkres(-obj.log_likelihood);
            
            obj.assignment_prob_ranks = [];
            obj.assignment_prob_probs = [];
            
            P = AutoId.sinkhorn_logspace(obj.log_likelihood, AutoId.iter_sinkhorn);
            
            assignment_det = zeros(1,size(P,1));
            for i=1:size(P,1)
                prow = [P(i,:) 1-nansum(P(i,:))];
                [sort_prob,ind_prob] = sort(-prow);
                ind_prob(ind_prob==length(prow))=-1;
                
                obj.assignment_prob_ranks(i,:) = ind_prob(1:obj.ncandidates);
                obj.assignment_prob_probs(i,:) = -sort_prob(1:obj.ncandidates);
                
                det = find(obj.assignments(i,:));
                
                if(isempty(det))
                    assignment_det(i)=-1;
                else
                    assignment_det(i)=det;
                end
            end
            
            for i=1:size(P,1)
                ind_det_in_probs = find(assignment_det(i)==obj.assignment_prob_ranks(i,:));
                rest_indices = setdiff(1:obj.ncandidates, ind_det_in_probs);
                rest_indices = rest_indices(1:obj.ncandidates-1);
                ranks = obj.assignment_prob_ranks(i,:);
                
                obj.assignment_prob_ranks(i, 1) = assignment_det(i);
                obj.assignment_prob_ranks(i, 2:obj.ncandidates)= ranks(rest_indices);
            end

            obj.assignments = obj.assignments;
            obj.assignment_prob_ranks = obj.assignment_prob_ranks;
            obj.assignment_prob_probs = obj.assignment_prob_probs;
            
        end
        
        function id(obj,im)
            % neurons that are already annotated
            annotations = im.get_annotations();
            annotation_confidences = im.get_annotation_confidences();
            [annotations{annotation_confidences<=0.5}] = deal('');
            
            annotated = [];
            annotated(:,1) = find(cellfun(@(x) ~isempty(x), annotations));
            annotated(:,2) = cellfun(@(x) find(strcmp(obj.atlas.(lower(im.bodypart)).N,x)), annotations(annotated(:,1)));
            
            
            % read features of the current image
            colors = im.get_colors_readout();
            colors = colors(:,[1 2 3]);
            
            % align the image to statistical atlas
            obj.align_to_atlas(colors, im.get_positions().*im.scale, ...
                                obj.atlas.(lower(im.bodypart)).model, annotated);
            
            obj.log_likelihood(annotated(:,1),:)= AutoId.min_log_likelihood;
            for i=1:size(annotated,1)
                obj.log_likelihood(annotated(i,1),annotated(i,2))= AutoId.max_log_likelihood;
            end
            
            % compute probabilistic/deterministic assignments
            obj.compute_assignments();
            
            % convert the assignments to names and update the information
            % in the image
            obj.add_to_image(im);
        end
        
        function update_id(obj, im, neuron_i, neuron)
            % update_auto_id changes the properties of auto_id given human has given 
            % the identity neuron to the neuron_i-th neuron (with the mp order).
            %it also updates the image object accordingly.
            obj.log_likelihood(neuron_i,:)= AutoId.min_log_likelihood;
            
            neuron_names = obj.atlas.(lower(im.bodypart)).N;
            
            if ~strcmp(neuron, 'NaN')
                obj.log_likelihood(neuron_i, strcmp(neuron_names, neuron))  = AutoId.max_log_likelihood;
            end
            
            obj.compute_assignments();
            
            % convert the assignments to names and update the information
            % in the image
            obj.add_to_image(im);
        end
        
        
        function align_to_atlas(obj, col, pos, model, annotated)
        % Automatic labeling of previously identified neurons based on mixture of
        % Gaussian model.
        %
        % Erdem & Amin
            h = waitbar(0,'Initialize ...');

            mu = nanmean(model.M,3); % should be replaced with model.mu;
            
            % initialization
            sigma_pooled_l = nancov(reshape(permute(model.M,[1 3 2]),[size(model.M,1)*size(model.M,3) size(model.M,2)]));
            sigma_pooled_u = nancov(reshape(permute(model.M0,[1 3 2]),[size(model.M0,1)*size(model.M0,3) size(model.M0,2)]));

            sigma = AutoId.lambda*model.Sigma+(1-AutoId.lambda)*sigma_pooled_l;
            
            sigma([1,2,3],[4,5,6]) = 0;
            sigma([4,5,6],[1,2,3]) = 0;
            
            model.Sigma = sigma;
            
            weights = ones(1,size(model.M,1));
            weights(annotated(:,2)) = 1/AutoId.annotation_weight;
%             Sigma_weighted = bsxfun(@times, repmat(eye(6),1,1,size(model.M,1)), reshape(weights,1,1,size(model.M,1)));
            
            Sigma_weighted = bsxfun(@times, sigma, reshape(weights,1,1,size(model.M,1)));
            
            init        = model;
            init.M      = model.M0;
            init.Sigma  = repmat(sigma_pooled_u,1,1,size(model.Sigma,3));
            
            P = AutoId.update_permutation(col, pos, init, annotated);
            
            beta_pos = AutoId.MCR_solver(mu(:,1:3), P'*[pos ones(size(pos,1),1)], Sigma_weighted(1:3,1:3,:));
            beta_col = AutoId.MCR_solver(mu(:,4:end), P'*[col ones(size(col,1),1)], Sigma_weighted(4:end,4:end,:));
            
            pos = [pos ones(size(pos,1),1)]*beta_pos;
            col = [col ones(size(col,1),1)]*beta_col;
                
            for iteration = 1: AutoId.iter_rwc
                try
                    waitbar(iteration/AutoId.iter_rwc,h,['Iteration ' num2str(iteration) '/' num2str(AutoId.iter_rwc)]);
                catch
                    break;
                end
                Z = arrayfun(@(i) mvnrnd(mu(i,:), model.Sigma(:,:,i)), 1:size(mu,1), 'UniformOutput', false);
                Z = vertcat(Z{:});

                P = AutoId.update_permutation(col,pos,model,annotated);
                
                beta_pos = AutoId.MCR_solver(Z(:,1:3), P'*[pos ones(size(pos,1),1)], Sigma_weighted(1:3,1:3,:));
                beta_col = AutoId.MCR_solver(Z(:,4:end), P'*[col ones(size(col,1),1)], Sigma_weighted(4:end,4:end,:));
                
                pos = [pos ones(size(pos,1),1)]*beta_pos;
                col = [col ones(size(col,1),1)]*beta_col;
            end
            
            obj.log_likelihood = - AutoId.pdist2_maha([pos col], mu, model.Sigma);
            
            try
                close(h);
            catch
                warning('ID is canceled.');
            end
        end
    end
end

