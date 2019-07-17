classdef AutoId < handle
    %AUTO_ID Auto.
    %
    %   An auto_id object is has a set of possible names (different if tail/head), the method used, a likelihood
    %   matrix, and a set of ids for each of observed neurons. 
    properties(Constant)
        % constants
        ncandidates         = 7;
        annotation_weight   = 0.5;
        
        iter_sinkhorn       = 2000;
        iter_rwc            = 10;
        
        lambda              = 0.5;
        
        min_log_likelihood  = -1e9
        max_log_likelihood  = 1e5
        
        theta               = 0: 0.25: 2*pi
    end
    
    properties
        
        atlas = getfield(load('atlas.mat'), 'atlas') % atlas data structure containing neuron names, positions, colors, and their covariances
        
        aligned % atlas aligned positions and colors
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
        
        function rotmat = rotmat(theta)
            rotmat = [cos(theta) -sin(theta);...
                      sin(theta) cos(theta)];
        end
        
        function pos = major_axis_align(pos,sign)
            [coeff,~,~,~,~,mu] = pca(pos);
            coeff(:,1:2) = coeff(:,1:2)*sign;
            pos = (pos-mu)*coeff;
        end
      
        function D = pdist2_maha(X,Y,Sigma)
            % mahalanobis distance between a set of points and the
            % components of a GMM
            D=zeros(size(X,1),size(Y,1));
            for i=1:size(Y,1)
                D(:,i)=pdist2(X,Y(i,:),'mahalanobis',Sigma(:,:,i));
            end
        end
        
        function vec_x = vec(x)
            % vectorize multi-dimensional array
            vec_x = x(:);
        end
        
        function beta = MCR_solver(Y,X,sigma)
            % MCR - Multiple covariance regression solver
            % Y - Target (n x d)
            % X - Source (n x p)
            % Sigma - covariances for each row of Y (d x d x n)
            % Solving: \sum_i \|Y_i - X_i \beta\|_{sigma_i}^2
            
            % variables
            d1 = size(Y,2); d2 = size(X,2); n = size(Y,1);
            
            % memory allocation
            A = zeros(d1*d2,d1*d2,n); B = zeros(d2,d1,n);
            
            % to learn about math look at the paper
            for i=1:size(Y,1)
                A(:,:,i) = kron(inv(sigma(:,:,i)),X(i,:)'*X(i,:));
                B(:,:,i) = X(i,:)'*Y(i,:)/sigma(:,:,i);
            end
            beta = reshape(nansum(A,3)\AutoId.vec(nansum(B,3)),[size(X,2),size(Y,2)]);
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
            M = model.mu;
            X = [positions colors];
            X(any(isnan(X),2),:) = [];
            
            log_likelihoodhat=zeros(size(M,1),size(X,1));
            order = randperm(size(M,3));
            for jj=1:length(order)
                j=order(jj);
                Y=M(:,:,j);
                idx=find(~any(isnan(Y),2));
                Y(any(isnan(Y),2),:)=[];
                D = AutoId.pdist2_maha(X,Y,model.sigma);
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
        
        function [col,pos,cost] = local_alignment(col,pos,model,theta,sgn,annotated)
            % local alignment between atlas and rotated points for given
            % theta
            POS = pos;
            
            % GMM for computing the likelihood of current local alignment
            gm = gmdistribution(model.mu,model.sigma);
            
            % standardize positions by axis aligning it
            pos = AutoId.major_axis_align(pos,sgn);
            
            % compute the rough transformation between atlas centers and
            % its axis-aligned version 
            [coeff_mu,~,~,~,~,mu_mu] = pca(model.mu(:,[1 2 3]));
            coeff_mu = coeff_mu*det(coeff_mu);
            
            % transform rotated positions based on the atlas transformation
            % to roughly align it
            pos(:,[2 3]) = pos(:,[2 3])*AutoId.rotmat(theta);
            pos = pos*coeff_mu'+mu_mu;
            
            % upweighting already annotated neurons
            weights = ones(1,size(model.mu,1));
            weights(annotated(:,2)) = 1/(size(model.mu,1)*AutoId.annotation_weight);
            
            sigma_weighted = bsxfun(@times, model.sigma, reshape(weights,1,1,size(model.mu,1)));
            
            for iteration = 1: AutoId.iter_rwc
                % sample from GMM components
                Z = arrayfun(@(i) mvnrnd(model.mu(i,:), model.sigma(:,:,i)), 1:size(model.mu,1), 'UniformOutput', false);
                Z = vertcat(Z{:});
                
                % update permutation
                P = AutoId.update_permutation(col,pos,model,annotated);
                
                % align point to sample
                beta_pos = AutoId.MCR_solver(Z(:,1:3), P'*[pos ones(size(pos,1),1)], sigma_weighted(1:3,1:3,:));
                beta_col = AutoId.MCR_solver(Z(:,4:end), P'*[col ones(size(col,1),1)], sigma_weighted(4:end,4:end,:));

                beta = AutoId.MCR_solver([pos(:,[1 2 3]) ones(size(pos,1),1)]*beta_pos,[POS ones(size(pos,1),1)],repmat(eye(3),1,1,size(pos,1)));
                
                % making sure that there are no mirror flips
                det_cost = sign(det(beta(1:3,1:3)));
                beta_pos = beta_pos*det_cost;
                
                % update positions and colors
                pos = [pos ones(size(pos,1),1)]*beta_pos;
                col = [col ones(size(col,1),1)]*beta_col;
            end
            
            % compute the cost based on log likelihood of current local
            % alignment and GMM distribution
            cost = nansum(log(gm.pdf([pos col])));
        end
        
    end
        
    
    methods
        function obj = AutoId()
            % load the statistical atlas
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
        
        function id(obj,im,varargin)
            if any(strcmp(varargin, 'atlas'))
                obj.atlas = varargin{find(strcmp(varargin, 'atlas'))+1};
            end
            
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
            obj.global_alignment(colors, im.get_positions().*im.scale, ...
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
            % it also updates the image object accordingly.
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
        
        
        
        function global_alignment(obj, col, pos, model, annotated)
        % Global alignment based on search in the theta space
        %
        % Amin & Erdem
            h = waitbar(0,'Initialize ...');
            
            cost = nan(2*length(AutoId.theta),1);
            colors = cell(2*length(AutoId.theta),1);
            positions = cell(2*length(AutoId.theta),1);
            
            for idx = 1:length(AutoId.theta)
                f(idx) = parfeval(@AutoId.local_alignment, 3, col, pos, model, AutoId.theta(idx), +1, annotated);
                f(length(AutoId.theta)+...
                  idx) = parfeval(@AutoId.local_alignment, 3, col, pos, model, AutoId.theta(idx), -1, annotated);
            end
            
            for idx=1:2*length(AutoId.theta)
                [job_idx,c,p,cc] = fetchNext(f);
                
                positions{job_idx}  = p;
                colors{job_idx}     = c;
                cost(job_idx)       = cc;
                try
                    waitbar(idx/(2*length(AutoId.theta)),h,[num2str(round((100.0*idx/(2*length(AutoId.theta))))),'%']);
                catch
                    break;
                end
            end
            
            try
                close(h);
            catch
                warning('Auto ID is canceled.');
            end
            
            [theta_idx] = find(cost==max(cost(:)));
            
            pos = positions{theta_idx};
            col = colors{theta_idx};

            obj.log_likelihood = -AutoId.pdist2_maha([pos col], model.mu, model.sigma);
            obj.aligned = [pos col];
        end
        
        function visualize(obj,im,varargin)
            if ~any(strcmp(varargin,'ax'))
                clf; hold on;
                set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
            end
            
            % original data
            model = obj.atlas.(lower(im.bodypart)).model;
            col = im.get_colors_readout(); col = col(:,[1 2 3]);
            pos = im.get_positions().*im.scale;
            
            % find transformation between original and aligned data
            beta = linsolve([obj.aligned ones(size(pos,1),1)],[pos col ones(size(pos,1),1)]);
            
            % move atlas based on the transformation
            model.mu = [model.mu ones(size(model.mu,1),1)]*beta;
            for i=1:size(model.sigma,3)
                model.sigma([1 2 3],[1 2 3],i) = beta([1 2 3],[1 2 3])'*model.sigma([1 2 3],[1 2 3],i)* beta([1 2 3],[1 2 3]);
            end
            
            % plot the result
            atlas_color = model.mu(:, [4 5 6]); atlas_color = (atlas_color)./(max(atlas_color));

            sizes = arrayfun(@(i) mean(eig(model.sigma([1,2,3],[1,2,3],i))), 1:size(model.sigma,3));
            sizes = 1000*sizes/max(sizes);
            
            if ~any(strcmp(varargin,'ax'))
                scatter3(model.mu(:,1),model.mu(:,2),model.mu(:,3), sizes, atlas_color, 'o', 'LineWidth', 3);
                text(model.mu(:,1),model.mu(:,2),model.mu(:,3), obj.atlas.(lower(im.bodypart)).N, 'Color', 'r');
                scatter3(pos(:,1),pos(:,2),pos(:,3), 100, col./max(col), '.');

                daspect([1,1,1]); grid on; set(gca,'color',[0.3,0.3,0.3]); view(-30,-20);
                drawnow;
            else
                ax = varargin{find(strcmp(varargin,'ax'))+1};
                z = varargin{find(strcmp(varargin,'z'))+1};
                
                mu_pixels = model.mu(:,1:3)./im.scale;
                current_z_indices = mu_pixels(:,3)>z-1.5 & mu_pixels(:,3)<z+1.5;
                scatter(ax,mu_pixels(current_z_indices,2),mu_pixels(current_z_indices,1),...
                    sizes(current_z_indices), atlas_color(current_z_indices,:), 'o', 'LineWidth', 3);
                text(ax,mu_pixels(current_z_indices,2),mu_pixels(current_z_indices,1),obj.atlas.(lower(im.bodypart)).N(current_z_indices), 'Color', 'r');
            end
        end

    end
    
    
end

