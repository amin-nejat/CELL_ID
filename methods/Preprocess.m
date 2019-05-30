classdef Preprocess < handle
    %PREPROCESS preprocesssing tools for NeuroPAL images
    
    properties
        
    end
    
    methods(Static)
        
        function obj = instance()
             persistent instance
             if isempty(instance)
                obj = AutoDetect();
                instance = obj;
             else
                obj = instance;
             end
        end
      
        function thresh = threshold_frame(video, method, sensitivity)
        % Thresholding multi-color volume
        %
        % Amin Nejat
            thresh = 0*video;

            for ch = 1: size(video, 4)
                vol = squeeze(video(:,:,:,ch,:));
                if strcmp(method, 'adaptive')
                    BW = imbinarize(vol/max(vol(:)), 'adaptive', 'Sensitivity', sensitivity);
                    vol(~BW) = 0;
                    vol = imclearborder(vol);
                    thresh(:,:,:,ch,:) = vol;
                elseif strcmp(method, 'global')
                    BW = imbinarize(vol/max(vol(:)), 'global');
                    vol(~BW) = 0;
                    vol = imclearborder(vol);
                    thresh(:,:,:,ch,:) = vol;
                elseif strcmp(method, 'otsu')
                    vol(vol < graythresh(vol(:))) = 0;
                    thresh(:,:,:,ch,:) = vol;
                elseif strcmp(method, 'fixed')
                    vol(vol < sensitivity) = 0;
                    thresh(:,:,:,ch,:) = vol;
                end
            end

        end
        
       function decimated = decimate_frame(video, new_size, method)
        % decimating a multi-dimensional array
        % DECIMATED = DECIMATE_FRAME(VIDEO, NEW_SIZE, METHOD) decimates the
        % multi-dimensional array VIDEO to size specified by NEW_SIZE. METHOD
        % options are:
        %   'linear'
        %   'nearest'
        %
        % See also RUN_FRAMES
        %
        % Amin Nejat
            new_size = round(new_size);
            new_size = new_size(1: 3);
            decimated = nan(new_size(1), new_size(2), new_size(3), size(video, 4));
            for ch = 1: size(video, 4)
                decimated(:, :, :, ch) = imresize3(video(:, :, :, ch), new_size, method);
            end
       end
       
       function hmvideo = histmatch_frame(video)
        % matching the histogram of the intensities for different colors to
        % normalize different color channels for better visualization
        % Amin Nejat
            hmvideo = zeros(size(video));
            exprand = exprnd(1, 1, numel(find(video > 0)));

            for ch = 1: size(video, 4)
                vol = squeeze(video(:,:,:,ch,:));
                posvol = vol(vol > 0);

                result = vol;
                result(vol > 0) = imhistmatchn(posvol(:)/max(posvol(:)), exprand/max(exprand(:)), 100000);
                hmvideo(:,:,:,ch,:) = reshape(result, size(vol));
            end
       end
        
       function volume = area_filter(volume, point, threshold, dim)
        % Heuristic filter for removing objects in a multi-color volume.
        %
        % Amin Nejat
            point = round(point);

            g = volume(:,:,:,dim);
            g(g < threshold) = 0;
            a = bwlabeln(g,26);

            mask = logical(0*a);
            mask(a == a(point(1), point(2), point(3))) = 1;
            mask = repmat(mask, 1,1,1,size(volume,4));
            volume(mask) = 0;

       end

        function zvideo = zscore_frame(video)
        % Z-scoring different color channels of the multi-color multi-dimensional
        % video. The colors must be in the 4th dimension. 
        %
        % Amin Nejat
            zvideo = zeros(size(video));
            for ch = 1: size(video, 4)
                data = double(video(:, :, :, ch, :));
                zvideo(:, :, :, ch, :) = reshape((data(:)-mean(data(:)))/std(data(:)), size(data));
            end
        end
        
        function rho = filter_frame(frame, filter)
        % applying a filter to multi-color RGB image
        % Amin Nejat
            rho = 0*frame;
            for ch = 1: size(frame, 4)
                D = squeeze(frame(:, :, :, ch));
                rho(:, :, :, ch) =  imfilter(D, filter, 'same');
            end
        end


    end
end

