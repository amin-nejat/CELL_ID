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
    methods
        
        
        function obj = AutoId(image, type)
            %class constructor, defines main properties of Auto_Id
            [neuron_names, LL, id_method] = auto_id(image, type); %call the main engine of auto_id
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
            %% add_to_image  function simply uptdates some properties of the image
            % structure that depend on auto_id

            image.add_deterministic_ids(obj.ids); %
            image.add_probabilistic_ids(obj.ids_prob);
            image.add_probabilistic_probs(obj.assignment_prob_probs)

            [~,rank_uncertain]=sort(obj.assignment_prob_probs(:,1));
            rank_uncertain_inv(rank_uncertain) = 1:length(rank_uncertain);
            image.add_ranks(rank_uncertain_inv');

        end
        
        function find_ids(obj)
            %% find_ids computes some secondary properties of auto_id object
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
            
           %% compute_assignments computes some secondary properties of auto_id object
            
            LL = obj.LL;
            
            [assignments,~]=munkres(-LL);
            
            ncandidates = 7;
            niter=2000;
            P = sinkhorn_logspace(LL, niter);
            
            
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
                %% update_auto_id_by_human changes the properties of auto_id given human has given 
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

