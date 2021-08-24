classdef NNDetect < handle
    %NNDETECT Detection of neural centers using neural nets
    
    properties(Constant)
        NN_model = 'NN_model.mat';
        stride = 128;
        crop_size = 128;
        pst_shape = [384,1024,48,4];
    end
    
    properties
        model_version
        model
    end
    
    methods
        function obj = NNDetect()
            %NNDETECT Construct an instance of this class
            %   Outputs:
            %       obj : an instance of this class
            
            % Load Neural Net Model.
            NN = load(Methods.NNDetect.NN_model);
            obj.model_version = NN.version;
            obj.model = NN.model;
        end
        
        function [pred_p] = predict_nn(obj,patches,pst_shape, titlestr)
            import Methods.*;
            
            % Setup the progress bar.
            wait_title = 'Detecting Neurons';
            wb = waitbar(0, {titlestr, 'Initializing ...'}, 'Name', wait_title);
            wb.Children.Title.Interpreter = 'none';
            
            num_patches = length(patches);
            pred = cell(num_patches,1);
            for i=1:num_patches
                
                % Update the progress bar.
                try
                    waitbar(i/num_patches, wb,...
                        {titlestr, ...
                        sprintf('%d%% completed ...', int16(100*i/num_patches))}, ...
                        'Name', wait_title);
                catch
                    break;
                end
                
                % Neural prediction.
                p = patches(i);
                pred{i} = predict(obj.model,p{1});
            end
            
            % Done.
            try
                close(wb);
            catch
                warning('The detection was canceled.');
                
                % Amin, can we use any predictions if the user cancels early?
                pred_p = [];
                return;
            end
            
            pred_ = zeros([NNDetect.pst_shape(1:3),1]);
            shape_ = floor(NNDetect.pst_shape(1:2)/NNDetect.stride);

            for i=1:shape_(1)
                for j=1:shape_(2)
                    pred_(NNDetect.stride*(i-1)+1:(i-1)*NNDetect.stride+NNDetect.crop_size,...
                        NNDetect.stride*(j-1)+1:(j-1)*NNDetect.stride+NNDetect.crop_size,...
                        :,:) = pred{(i-1)*shape_(2)+j};
                end
            end
            
            pred_p = NNDetect.pad_to_shape(pred_,pst_shape(1:3));
        end
    end
    
    methods(Static)
        
        function obj = instance()
            persistent instance
            if isempty(instance)
                obj = Methods.NNDetect();
                instance = obj;
            else
                obj = instance;
            end
        end
        
        function [supervoxels, params] = detect(titlestr, data)
            %DETECT Use a neural network to detect neurons in the image volume.
            % Inputs:
            %   titlestr = a title string for the progress bar
            %   data   = the image volume to search for neurons
            % Outputs:
            %   supervoxels = struct containing the locations (positions), shape
            %    parameters (covariances), colors (colors), basline color
            %    (baseline), truncation values (trunc), and readout color from
            %    the image (color_readout)
            %   params      = struct contining the parameters chosen for running the
            %    algorithm; parameters are half size of a neuron in microns
            %    (hnsz), number of objects (k), threshold for minimum eigenvalue
            %    (min_eig_thresh), exclusion radius (exclusion_radius)
            import Methods.NNDetect;
            
            % Detect the neurons.
            data_p = NNDetect.pad_to_shape(data, NNDetect.pst_shape);
            patches = NNDetect.create_patches(data_p);
            pred_p = NNDetect.instance().predict_nn(patches,size(data), titlestr);
            
            % Do we have any neuron predictions?
            if isempty(pred_p)
                supervoxels = [];
                params = [];
                return;
            end
            
            % Store the detected neurons.
            pred_p = (pred_p > 0.5)*1;
            CC = bwconncomp(pred_p,26);
            stats = regionprops3(CC,'Centroid');
            centers = stats.Centroid;
            centers = centers(:,[2,1,3]);
            
            % Store the detected neuron centers.
            supervoxels = [];
            supervoxels.positions = centers;
            
            % Store the detected neuron covariances.
            supervoxels.covariances = zeros(size(centers,1),3,3);
            supervoxels.covariances(:,1,1) = 10;
            supervoxels.covariances(:,2,2) = 10;
            supervoxels.covariances(:,3,3) = 10;
            
            % Store the detected neuron colors.
            supervoxels.color = zeros(size(centers,1),4);
            supervoxels.color_readout = zeros(size(centers,1),4);

            for n=1:size(centers,1)
                pos = round(centers(n,:));
                supervoxels.color(n,:) = data(pos(1),pos(2),pos(3),:);
                supervoxels.color_readout(n,:) = data(pos(1),pos(2),pos(3),:);
            end
            
            supervoxels.baseline = zeros(size(centers,1),4);
            supervoxels.truncation = zeros(size(centers,1),1);
            
            % Store the neural detection parameters.
            params                  = [];
            params.k                = size(supervoxels.color,1);
            params.detect_scale     = 0;
            params.hnsz             = [10,10,3];
            params.min_eig_thresh   = 0; % microns
            params.exclusion_radius = 0; % microns
        end
        
        function BatchDetect(file, worm, num_neurons)
            %BATCHDETECT Batch detect neurons.
            
            % Setup the conversion progress bar.
            % Note: windows wants the interpreter off from the beginning.
            wait_title = 'Converting Image';
            wb = waitbar(0, 'Converting ...', 'Name', wait_title);
            wb.Children.Title.Interpreter = 'none';
            waitbar(0, wb, {file, 'Converting ...'}, 'Name', wait_title);
            
            % Open the image file.
            try
                [data, info, prefs, ~, mp, ~, np_file, id_file] = ...
                    DataHandling.NeuroPALImage.open(file);
            catch
                % For now, don't throw exceptions from threads.
                warning('Cannot read: "%s"', file);
                return;
            end
            
            % Done converting.
            try
                close(wb);
            catch
                warning('Image conversion was canceled.');
            end
            
            % Save the worm info.
            save(np_file, 'worm', '-append');
           
            % Setup the preprocessing progress bar.
            % Note: windows wants the interpreter off from the beginning.
            wait_title = 'Preprocessing Image';
            wb = waitbar(0, 'Preprocessing ...', 'Name', wait_title);
            wb.Children.Title.Interpreter = 'none';
            waitbar(0, wb, {file, 'Preprocessing ...'}, 'Name', wait_title);
            
            % Preprocess the colors.
            data_RGBW = double(data(:,:,:,prefs.RGBW(~isnan(prefs.RGBW))));
            data_zscored_raw = Methods.Preprocess.zscore_frame(double(data_RGBW));
            
            % Remove artifacts.
            [~, mask] = Methods.Preprocess.filter_gut_lysosomes(data_zscored_raw);
            mask = repmat(mask,1,1,1,length(prefs.RGBW(~isnan(prefs.RGBW))));
            data_zscored = data_zscored_raw;
            data_zscored(mask) = 0;
            
            % Done preprocessing.
            try
                close(wb);
            catch
                warning('Image preprocessing was canceled.');
            end
            
            % Detect the neurons.
            mp.k = num_neurons;
            [sp, mp] = Methods.NNDetect.instance().detect(file, data_zscored);
            
            % Save the neurons.
            version = Program.ProgramInfo.version;
            mp_params = mp;
            neurons = Neurons.Image(sp, worm.body, 'scale', info.scale');
            Methods.Utils.removeNearbyNeurons(neurons, 2, 2);
            save(id_file, 'version', 'neurons', 'mp_params');
        end
        
        function [data_zscored,patches] = prepare_image(file)
            import Methods.*;
            
            [data, ~, prefs, ~, ~, ~, ~, ~] = DataHandling.NeuroPALImage.open(file);
            data_zscored = Methods.Preprocess.zscore_frame(data);
            data_zscored = data_zscored(:,:,:,prefs.RGBW(1:4));
            data_p = NNDetect.pad_to_shape(data_zscored,NNDetect.pst_shape);
            patches = NNDetect.create_patches(data_p);
        end
        
        function [data_p] = pad_to_shape(data,pst_shape)
            pre_shape = size(data);
            pds = ceil((pre_shape-pst_shape)/2);
            pds(pds<0) = 0;
            
            sls = ceil(-(pre_shape-pst_shape)/2);
            sls(sls < 0) = 0;
            data_p = padarray(data,sls,'symmetric');
            
            data_p = data_p(pds(1)+1:pds(1)+pst_shape(1),...
                pds(2)+1:pds(2)+pst_shape(2),...
                pds(3)+1:pds(3)+pst_shape(3),:);
            
        end
        
        function patches = create_patches(data)
            import Methods.*;
            
            patches = {};
            
            shape_ = size(data);
            
            for i=1:floor(shape_(1)/NNDetect.stride)
                if (i-1)*NNDetect.stride+NNDetect.crop_size <= shape_(1)
                    for j=1:floor(shape_(2)/NNDetect.stride)
                        if (j-1)*NNDetect.stride+NNDetect.crop_size <= shape_(2)
                            patches{end+1} =  data(NNDetect.stride*(i-1)+1:(i-1)*NNDetect.stride+NNDetect.crop_size, ...
                                NNDetect.stride*(j-1)+1:(j-1)*NNDetect.stride+NNDetect.crop_size,:,:);
                            
                        end
                    end
                end
            end
        end
        
        function visualize(data,centers,pred_p,save,varargin)
            figure('Renderer', 'painters', 'Position', [10 10 900 600]);
            
            subplot(211);
            image(double(max(pred_p,[],3))*256);
            set(gca,'XTickLabel',[]);
            set(gca,'YTickLabel',[]);
            
            subplot(212); hold on;
            image(squeeze(max(data(:,:,:,1:3),[],3))/20);
            scatter(centers(:,1),centers(:,2),10,'MarkerFaceColor','k','MarkerEdgeColor','w',...
                'LineWidth',1);
            xlim([1,size(data,2)]);
            ylim([1,size(data,1)]);
            set(gca,'XTickLabel',[]);
            set(gca,'YTickLabel',[]);
            set(gca,'YDir','reverse');
            
            if save
                saveas(gcf,[varargin{1},'.png']);
                saveas(gcf,[varargin{1},'.pdf']);
                close('all');
            end
        end
        
        function [centers] = post_process(pred_p)
            pred_p = (pred_p > 0.5)*1;
            
            CC = bwconncomp(pred_p,26);
            stats = regionprops3(CC,'Centroid');
            centers = stats.Centroid;
        end
    end
end

