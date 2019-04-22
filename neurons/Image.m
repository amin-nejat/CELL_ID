classdef Image < handle
    %Image A list of neurons in a certain body part
    %   Detailed explanation goes here

    properties
        neurons = Neuron.empty; % a list of the instances of Neuron class -> see Neuron.m
        bodypart % a string consisting the name of the worm's body part
        meta_data = containers.Map(); % key, value pairs for intermediate analysis
        scale = ones(1,3); % (x,y,z) scale
    end

    methods
        function obj = Image(superpixels, bodypart, varargin)
            %Image Construct an instance of this class.
            %   superpixel: Matlab struct superpixels with variables mean, cov,
            %   color, basline, and potentially ids, rank, probabilistic
            %   ids, probabilistic probs, ... .
            %   bodypart: A string that represents which body part the
            %   current instance of this class corresponds to, examples are
            %   'head' and 'tail'.
            %   [scale]: optional image scale (x,y,z).

            % No neurons.
            obj.bodypart = bodypart;
            if isempty(superpixels)
                obj.neurons = [];
                return;
            end

            % Create the neurons.
            for i=1:length(superpixels.mean)
                obj.neurons(i) = Neuron(sub_sp(superpixels,i));
            end

            % Setup the image scale.
            if ~isempty(varargin)
                obj.scale = varargin{1};
            end
        end

        function add_neuron(obj, volume, position, nsz, trunc)
            %ADD_NEURON Adds a neuron to the list of neurons in obj.neurons
            %   by running one iteration of Matching Pursuit.
            %   volume: the full z-scored image.
            %   position: the location in the neighbourhood of which the
            %   neuron should be added.
            %   nsz: the window size around the patch for the  new neuron.
            %   trunc: the truncation value of the Gaussian function used
            %   for fitting.
            color = squeeze(volume(round(position(1)),round(position(2)),round(position(3)),:))';
            bpatch = subcube(volume, round(position), nsz);
            if isKey(obj.meta_data, 'auto_detect') && obj.meta_data('auto_detect')
                [~, sp, ~] = fit_gaussian(double(bpatch), size(volume), color, nsz, trunc, position);
            else
                sp = [];
                sp.mean = position;
                sp.color = color;
                sp.baseline = [0,0,0,0];
                sp.cov = diag(nsz);
            end

            if isempty(obj.neurons)
                obj.neurons = Neuron(sp);
            else
                obj.neurons(end+1) = Neuron(sp);
            end
        end

        function rot_image = rotate_X_180(obj, image)
            %ROTATE_X_180 Rotate the image 180 degrees around the x-axis.
            %   image: the image to rotate
            %   rot_image: the image rotated 180 degrees around the x-axis.
            rot_image = image(:,end:-1:1,end:-1:1,:,:);
            for i=1:length(obj.neurons)
                position =  obj.neurons(i).position;
                obj.neurons(i).position(2) = size(image,2) - position(2) + 1;
                obj.neurons(i).position(3) = size(image,3) - position(3) + 1;
            end
        end

        function rot_image = rotate_Y_180(obj, image)
            %ROTATE_Y_180 Rotate the image 180 degrees around the y-axis.
            %   image: the image to rotate
            %   rot_image: the image rotated 180 degrees around the y-axis.
            rot_image = image(end:-1:1,:,end:-1:1,:,:);
            for i=1:length(obj.neurons)
                position =  obj.neurons(i).position;
                obj.neurons(i).position(1) = size(image,1) - position(1) + 1;
                obj.neurons(i).position(3) = size(image,3) - position(3) + 1;
            end
        end

        function rot_image = rotate(obj, image, rot)
            %ROT_IMAGE Rotates the image using parameters specified in rot
            %   image: the full image that needs to be rotated
            %   rot: a 4x4 rotation matrix or a 1x3 vector consisting
            %   rotation angles. The rotation will be applied to the full
            %   image and the neurons one by one.
            %   rot_image: the rotated image.
            if isvector(rot)
                rot_mat = makehgtform('xrotate',rot(1),'yrotate',rot(2),'zrotate',rot(3));
            else
                rot_mat = rot;
            end

            rot_image = [];
            sz = size(image);
            Rin = imref3d(sz(1:3));

            Rin.XWorldLimits = Rin.XWorldLimits-2*mean(Rin.XWorldLimits);
            Rin.YWorldLimits = Rin.YWorldLimits-2*mean(Rin.YWorldLimits);
            Rin.ZWorldLimits = Rin.ZWorldLimits-2*mean(Rin.ZWorldLimits);

            for c=1:size(image,4)
                rot_image(:,:,:,c) = imwarp(image(:,:,:,c,:),Rin,affine3d(rot_mat),'cubic');
            end

            newsz = size(rot_image);

            for i=1:length(obj.neurons)
                obj.neurons(i).rotate(rot_mat,sz([1,2,3]),newsz([1,2,3]))
            end
        end

        function num = num_neurons(obj)
            %NUM_NEURONS the number of neurons in the image
            num = length(obj.neurons);
        end

        function num = num_user_id_neurons(obj)
            %NUM_USER_ID_NEURONS the number of user ID'd neurons in the image
            num = sum(arrayfun(@(x) ~isempty(x.annotation), obj.neurons));
        end

        function is_fully_annotated = is_all_annotated(obj)
            %IS_ALL_ANNOTATED do all neurons have user annotations?
            is_fully_annotated = obj.num_neurons() == obj.num_user_id_neurons();
        end

        function names = user_id_neuron_names(obj)
            %USER_ID_NEURON_NAMES get a list of the user ID neuron names
            names = {};
            for i = 1:length(obj.neurons)
                if ~isempty(obj.neurons(i).annotation)
                    names{end+1} = obj.neurons(i).annotation;
                end
            end
        end

        function neurons = user_id_neurons(obj)
            %USER_ID_NEURONS get a list of the user ID neurons
            id_neurons = arrayfun(@(x) ~isempty(x.annotation), obj.neurons);
            neurons = obj.neurons(id_neurons);
        end

        function [neuron, i] = find_user_id_neuron(obj, name)
            %FIND_USER_ID find a neuron by user ID name
            neuron = [];
            for i = 1:length(obj.neurons)
                if strncmp(name, obj.neurons(i).annotation, length(name))
                    neuron = obj.neurons(i);
                    return;
                end
            end
            i = [];
        end

        function [neuron, i] = nearest_unannotated(obj, position)
            %NEAREST_UNANNOTATED_Z find the nearest unannoted neuron in z

            % Find the unannotated neurons.
            unIDd_i = find(arrayfun(@(x) isempty(x.annotation), obj.neurons));
            positions = round(vertcat(obj.neurons(unIDd_i).position));

            % No neurons found.
            neuron = [];
            i = [];
            if isempty(positions)
                return;
            end
            
            % Correct the positions for image scale.
            if size(position,1) > size(position,2)
                position = position';
            end
            position = position .* obj.scale';
            positions = positions .* obj.scale';
            
            % Find the nearest unannotated neuron in (x,y,z).
            [~,min_i] = min(sum((positions - position).^2,2));
            i = unIDd_i(min_i);
            neuron = obj.neurons(i);
        end

        function [neuron, i] = nearest_unannotated_z(obj, position)
            %NEAREST_UNANNOTATED_Z find the nearest unannoted neuron in z
            
            % Find the unannotated neurons.
            position = round(position);
            unIDd_i = find(arrayfun(@(x) isempty(x.annotation), obj.neurons));
            positions = round(vertcat(obj.neurons(unIDd_i).position));
            
            % No neurons found.
            neuron = [];
            i = [];
            if isempty(positions)
                return;
            end
            
            % Find the nearest unannotated neuron in z.
            [~, min_z_i] = min(abs(positions(:,3) - position(3)));
            min_z = positions(min_z_i,3);
            min_z_i = find(positions(:,3) == min_z);
            
            % Find the nearest unannotated neuron in (x,y).
            [~, min_xy_i] = min(sum(positions(min_z_i,1:2).^2) - sum(position(1:2).^2));
            i = unIDd_i(min_z_i(min_xy_i));
            neuron = obj.neurons(i);
        end
        
        function annotations = get_annotations(obj)
            %GET_ANNOTATIONS getter of neuron annotations.
            annotations = vertcat({obj.neurons.annotation});
        end

        function is_annotations_on = get_is_annotations_on(obj)
            %GET_IS_ANNOTATIONS_ON getter of neuron is_annotation_on s.
            is_annotations_on = vertcat(obj.neurons.is_annotation_on);
        end

        function annotation_confidences = get_annotation_confidences(obj)
            %GET_ANNOTATION_CONFIDENCES getter of neuron annotation_confidence s.
            annotation_confidences = vertcat(obj.neurons.annotation_confidence);
        end

        function add_deterministic_ids(obj, deterministic_ids)
            %ADD_DETERMINISTIC_IDS setter of neuron deterministic_id s.
            for i=1:length(obj.neurons)
                obj.neurons(i).deterministic_id = deterministic_ids{i};
            end
        end

        function deterministic_ids = get_deterministic_ids(obj)
            %GET_DETERMINISTIC_IDS getter of neuron deterministic_id s.
            deterministic_ids = vertcat({obj.neurons.deterministic_id});
        end

        function add_probabilistic_ids(obj, probabilistic_ids)
            %ADD_PROBABILISTIC_IDS stter of neuron probabilistic_id s.
            for i=1:length(obj.neurons)
                obj.neurons(i).probabilistic_ids = probabilistic_ids(i, :);
            end
        end

        function probabilistic_ids = get_probabilistic_ids(obj)
            %GET_PROBABILISTIC_IDS getter of neuron probabilistic_id s.
            probabilistic_ids = vertcat(obj.neurons.probabilistic_ids);
        end

        function add_probabilistic_probs(obj, probabilistic_probs)
            %ADD_PROBABILISTIC_PROBS setter of neuron probabilistic_prob s.
            for i=1:length(obj.neurons)
                obj.neurons(i).probabilistic_probs = probabilistic_probs(i, :);
            end
        end

        function probabilistic_probs = get_probabilistic_probs(obj)
            %GET_PROBABILISTIC_PROBS getter of neuron probabilistic_prob s.
            probabilistic_probs = vertcat(obj.neurons.probabilistic_probs);
        end

        function add_ranks(obj, ranks)
            %ADD_RANKS setter of neuron rank s.
            for i=1:length(obj.neurons)
                obj.neurons(i).rank = ranks(i);
            end
        end

        function ranks = get_ranks(obj)
            %GET_RANKS getter of neuron rank s.
            ranks = vertcat(obj.neurons.rank);
        end


        function colors = get_marker_colors(obj)
            %GET_MARKER_COLORS getter of neuron marker_color s.
            colors = zeros(length(obj.neurons), 3);
            for i=1:length(obj.neurons)
                colors(i,:) = obj.neurons(i).get_marker_color();
            end
        end

        function sizes = get_marker_sizes(obj)
            %GET_MARKER_SIZES getter of neuron marker_size s.
            sizes = zeros(length(obj.neurons), 1);
            for i=1:length(obj.neurons)
                sizes(i) = obj.neurons(i).get_marker_size();
            end
        end

        function sizes = get_line_sizes(obj)
            %GET_LINE_COLORS getter of neuron line_size s.
            sizes = zeros(length(obj.neurons), 1);
            for i=1:length(obj.neurons)
                sizes(i) = obj.neurons(i).get_line_size();
            end
        end

        function image = get_3d_shape(obj,sz,nsz,trunc)
            %GET_MARKER_COLORS getter of neuron 3d reconstruction s.
            %   sz: size of the full reconstructed image, which is the same
            %   size of the original image used for fitting.
            %   nsz: size of a neuron used for reconstruction.
            %   trunc: Gaussian truncation value used for reconstruction.
            %   image: full reconstruction of image using truncated
            %   Gaussian functions and stitching them together.
            image = zeros(sz);
            for i=1:length(obj.neurons)
                image = superpose(image, round(obj.neurons(i).position), obj.neurons(i).get_3d_reconstruction(nsz,trunc));
            end
        end

        function positions = get_positions(obj)
            %GET_POSITIONS getter of neuron position s.
            positions = vertcat(obj.neurons.position);
        end

        function colors = get_colors(obj)
            %GET_COLORS getter of neuron color s.
            colors = vertcat(obj.neurons.color);
        end

        function baselines = get_baselines(obj)
            %GET_BASELINES getter of neuron baseline s.
            baselines = vertcat(obj.neurons.baseline);
        end

        function covariances = get_covariances(obj)
            %GET_COVARIANCES getter of neuron covariance s.
            covariances = permute(cat(3, obj.neurons.covariance), [3,1,2]);
        end

        function add_meta_data(obj, key, value)
            %ADD_META_DATA adding a (key,value) pair in the meta_data data
            %structure for intermediate analysis, examples are LL (log
            %likelihood) etc.
            obj.meta_data(key) = value;
        end

        function value = get_meta_data(obj, key)
            %GET_META_DATA returns the value paired with key.
            value = obj.meta_data(key);
        end


        function sp = to_superpixel(obj)
            %TO_SUPERPIXEL coverts the neurons to a superpixel data
            %structure for data loading and storing and for interfacing
            %with different softwares and programming languages.
            sp = [];
            sp.mean = obj.get_positions();
            sp.color = obj.get_colors();
            sp.baseline = obj.get_baselines();
            sp.cov = obj.get_covariances();
            sp.probabilistic_ids = obj.get_probabilistic_ids();
            sp.deterministic_id = obj.get_deterministic_ids()';
            sp.annotation = obj.get_annotations()';
            sp.is_annotation_on = obj.get_is_annotations_on();
            sp.probabilistic_probs = obj.get_probabilistic_probs();
            sp.annotation_confidence = obj.get_annotation_confidences();
        end
    end
end
