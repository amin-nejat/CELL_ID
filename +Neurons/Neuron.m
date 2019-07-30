classdef Neuron < handle
    %NEURON A neuron.
    %
    %   A neuron has a position within an image, color, list of potential
    %   IDs, probabilities associated with these IDs, and methods to draw
    %   it within a GUI.

    properties (Access = private)
        brain
    end

    properties
        % Neuron position & color.
        position % neuron pixel position (x,y,z)
        color % neuron color based on fitting (R,G,B,W,...), W = white channel, values=[0-255]
        color_readout % neuron color based on readout from the image
        baseline % baseline noise values (R,G,B,W,...), values=[-1,1]
        covariance % 3x3 covariance matrix that's fit to the neuron
        truncation = 0 % Gaussian truncation value that defines the sharpness of the edge of neuron
        aligned_xyzRGB % the neuron's position & color, aligned to the global model
                
        % User neuron ID.
        annotation = '' % neuron user selected annotation
        is_annotation_on = NaN % is the neuron's annotation ON, OFF, or neither (empty)
        annotation_confidence = -1 % user confidence about annotation
        
        % Auto neuron ID.
        deterministic_id  % neuron ID assigned by the deterministic model
        probabilistic_ids % neuron IDs listed by descending probability
        probabilistic_probs % neuron ID probabilities
        rank % ranks of the neuron based on the confidence assigned by sinkhorn algorithm
        
        % GUI properties.
        % MOVE TO IMAGE!!!!
        is_selected = false % GUI related parameter specifying if the neuron is selected in the software or not
        MARKER_SIZE_NOT_SELECTED = 40 % unselected neuron marker size
        MARKER_SIZE_SELECTED = 200 % selected neuron marker size
        LINE_SIZE_NOT_SELECTED = 1 % unselected neuron line size
        LINE_SIZE_SELECTED = 4 % selected neuron line size
    end

    methods (Static)
        function neuron = unmarshall(sp, i)
            % UNMARSHALL Construct the neuron by unmarshalling it.
            % *** LEGACY FUNCTION. DEPRECATED!!!
            
            % Construct the neuron.
            neuron = Neurons.Neuron;
            
            % Neuron position & color.
            neuron.position = sp.positions(i,:);
            neuron.color = sp.color(i,:);
            neuron.color_readout = sp.color_readout(i,:);
            neuron.baseline = sp.baseline(i,:);
            neuron.covariance = sp.covariances(i,:,:);
            if isfield(sp, 'truncation')
                neuron.truncation = sp.truncation(i,:);
            end
            if isfield(sp, 'aligned_xyzRGB')
                if ~isempty(sp.aligned_xyzRGB) && ...
                        i <= size(sp.aligned_xyzRGB,1)
                    neuron.aligned_xyzRGB = sp.aligned_xyzRGB(i,:);
                end
            end
            
            % User neuron ID.
            if isfield(sp, 'annotation')
                neuron.annotation = sp.annotation{i};
                neuron.is_annotation_on = sp.is_annotation_on(i);
                neuron.annotation_confidence = sp.annotation_confidence(i);
            end
            
            % Auto neuron ID.
            if isfield(sp, 'deterministic_id')
                neuron.deterministic_id = sp.deterministic_id{i};
                if size(sp.probabilistic_ids,1) > 0 && ...
                        i <= size(sp.probabilistic_ids,1)
                    neuron.probabilistic_ids = sp.probabilistic_ids(i,:);
                end
                if size(sp.probabilistic_probs,1) > 0 && ...
                        i <= size(sp.probabilistic_probs,1)
                    neuron.probabilistic_probs = sp.probabilistic_probs(i,:);
                end
                if isfield(sp, 'rank')
                    if ~isempty(sp.rank) && i <= length(sp.rank)
                        neuron.rank = sp.rank(i);
                    end
                end
            end
        end
    end
    
    methods
        function annotate(obj, name, confidence, is_on)
            % ANNOTATE annotate the neuron.
            %   name: the neuron name
            %   confidence: the user confidence

            % Remove the user annotation.
            if isempty(name)
                obj.annotation = '';
                obj.annotation_confidence = -1;
                obj.is_annotation_on = nan;

            % Annotate the neuron.
            else
                obj.annotation = name;
                obj.annotation_confidence = confidence;
                obj.is_annotation_on = is_on;
            end
        end

        function delete_annotation(obj)
            %DELETE_ANNOTATION delete the user ID.
            obj.annotate('', 0, nan);
        end
        
        function delete_model_ID(obj)
            %DELETE_MODEL_ID delete the model-predicted ID.
            obj.deterministic_id = [];
            obj.probabilistic_ids = [];
            obj.probabilistic_probs = [];
            obj.rank = [];
            obj.aligned_xyzRGB = [];
        end
        
        function rotate(obj, rot, sz, newsz)
            %ROTATE rotates the neuron and translates it according to new
            %image space.
            %   rot: a 4x4 rotation matrix or a 1x3 vector consisting
            %   rotation angles.
            %   sz: size of the image before rotation used to translate
            %   before rotation.
            %   newsz: size of the image after rotation used to translate
            %   after rotation.

            if isvector(rot)
                rot_mat = makehgtform('xrotate',rot(1),'yrotate',rot(2),'zrotate',rot(3));
            else
                rot_mat = rot;
            end

            obj.position([2,1,3]) = (obj.position([2,1,3])-((sz([2,1,3])+1)/2))*rot_mat(1:3,1:3)+((newsz([2,1,3])+1)/2);
            obj.covariance([2,1,3],[2,1,3]) = rot_mat(1:3,1:3)*obj.covariance([2,1,3],[2,1,3])*rot_mat(1:3,1:3)';
        end

        function addToBrain(obj, b)
            %ADDTOBRAIN Add this neuron to a brain image.
            obj.brain = b;
        end

        function remove(obj)
            %REMOVE Remove this neuron.
            if ~isempty(obj.brain)
                obj.brain.remove(obj);
            end
        end

        function color = get_marker_color(obj)
            %GET_MARKER_COLORS specifies the marker color of the current
            %neuron according to the annotation confidence.
            if obj.is_selected
                color = [1,1,1]; % neuron selected (white)
            elseif obj.annotation_confidence == 0
                color = [1,0,0]; % user added neuron but no ID yet (red)
            elseif obj.annotation_confidence == 1
                color = [0,1,0]; % user ID'd neuron as 100% correct (green)
            elseif obj.annotation_confidence == 0.5
                color = [1,1,0]; % user ID'd neuron as low probability (yellow)
            elseif obj.annotation_confidence == -1
                color = [1,0.5,0]; % model ID'd neuron (orange)
            end
        end

        function size = get_marker_size(obj)
            %GET_MARKER_SIZE specifies the marker size of the current
            %neuron according to whether it is selected in the software or
            %not.
            if ~obj.is_selected
                size = obj.MARKER_SIZE_NOT_SELECTED;
            else
                size = obj.MARKER_SIZE_SELECTED;
            end
        end

        function size = get_line_size(obj)
            %GET_MARKER_SIZE specifies the line size of the current
            %neuron according to whether it is selected in the software or
            %not.
            if ~obj.is_selected
                size = obj.LINE_SIZE_NOT_SELECTED;
            else
                size = obj.LINE_SIZE_SELECTED;
            end
        end

        function recon = get_3d_reconstruction(obj,nsz)
            %GET_3D_RECONSTRUCTION reconstructs the 3D colored shape of the
            %current neuron according to its properties.
            %   nsz: size of the reconstructed image
            %   trunc: truncation value used for reconstruction.
            %   recon: (X,Y,Z,C) 4D double array.
            recon = Utils.simulate_gaussian(2*nsz+1, ... % size
                nsz+1+obj.position-round(obj.position), ... % center
                obj.covariance, ... % cov
                obj.color, ... % colors
                zeros(size(obj.color)), ... % baseline mean
                obj.truncation); % truncation
        end
    end
end
