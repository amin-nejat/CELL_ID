classdef Image < handle
    %Image A list of neurons in a certain body part
    %   Detailed explanation goes here
    
    properties
        neurons % a list of the instances of Neuron class -> see Neuron.m
        bodypart % a string consisting the name of the worm's body part 
        meta_data = containers.Map(); % key, value pairs for intermediate analysis
    end
    
    methods
        function obj = Image(superpixels, bodypart)
            %Image Construct an instance of this class.
            %   superpixel: Matlab struct superpixels with variables mean, cov,
            %   color, basline, and potentially ids, rank, probabilistic
            %   ids, probabilistic probs, ... .
            %   bodypart: A string that represents which body part the
            %   current instance of this class corresponds to, examples are
            %   'head' and 'tail'.
            
            if isempty(superpixels)
                obj.neurons = [];
                obj.bodypart = bodypart;
                return;
            end
            
            for i=1:length(superpixels.mean)
                neurons(i) = Neuron(sub_sp(superpixels,i));
            end
            obj.neurons = neurons;
            obj.bodypart = bodypart;
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
                sp.baseline = [0,0,0];
                sp.cov = diag(nsz);
            end
            
            if isempty(obj.neurons)
                obj.neurons = Neuron(sp);
            else
                obj.neurons(end+1) = Neuron(sp);
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
        
        function annotations = get_annotations(obj)
            %GET_ANNOTATIONS getter of neuron annotations.
            annotations = vertcat({obj.neurons.annotation});
        end
        
        function annotation_confidences = get_annotation_confidences(obj)
            %GET_ANNOTTION_CONFIDENCES getter of neuron annotation_confidence s.
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
            sp.probabilistic_probs = obj.get_probabilistic_probs();
            sp.annotation_confidence = obj.get_annotation_confidences();
        end
    end
end

