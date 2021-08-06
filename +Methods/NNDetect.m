classdef NNDetect < handle
    %NNDETECT Detection of neural centers using neural nets
    
    properties(Constant)
        stride = 128;
        crop_size = 128;
        pst_shape = [384,1024,48,4];
    end
    
    properties
        model
    end
    
    methods
        function obj = NNDetect(file)
            %NNDETECT Construct an instance of this class
            %   Inputs:
            %       file : file address of the pre-trained tensorflow model
            %
            %   Outputs:
            %       obj : an instance of this class
            
            obj.model = importKerasNetwork(file);
        end
        
        function [pred_p] = predict_nn(obj,patches,pst_shape, titlestr)
            import Methods.*;
            
            % Setup the progress bar.
            wait_title = 'Detecting Neurons';
            wb = waitbar(0, {titlestr, 'Initializing ...'}, 'Name', wait_title);
            wb.Children.Title.Interpreter = 'none';
        
            pred = {};
            num_patches = length(patches);
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
                pred{end+1} = predict(obj.model,p{1});
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

