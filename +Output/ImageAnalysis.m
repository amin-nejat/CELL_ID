classdef ImageAnalysis
    %IMAGEANALYSIS Methods to analyze NeuroPAL images.
    
    % Analysis constants.
    properties (Constant, Access = public)
        GFP_bg_size = 5; % GFP background size (pixels squared)
    end
    
   % Public methods.
    methods (Static)
        function saveID2CSV(csvfile, prefs, data, data_zscored, neurons, um_scale)
            %SAVEID2CSV save the IDs, position, & colors to a CSV file.

            % Measure the 8 image corners to determine an appropriate
            % GFP background threshold.
            GFP_bg_size = Output.ImageAnalysis.GFP_bg_size;
            GFP_i = prefs.GFP;
            GFP_image = double(squeeze(data(:,:,:,GFP_i)));
            GFP_image = (GFP_image - mean(GFP_image(:))) / std(GFP_image(:));
            [x,y,z] = size(GFP_image);
            x1 = 1:GFP_bg_size;
            x2 = (x - GFP_bg_size + 1):x;
            y1 = 1:GFP_bg_size;
            y2 = (y - GFP_bg_size + 1):y;
            z1 = 1;
            z2 = z;
            corner(1) = median(GFP_image(x1,y1,z1), 'all');
            corner(2) = median(GFP_image(x2,y1,z1), 'all');
            corner(3) = median(GFP_image(x1,y2,z1), 'all');
            corner(4) = median(GFP_image(x1,y1,z2), 'all');
            corner(5) = median(GFP_image(x2,y2,z1), 'all');
            corner(6) = median(GFP_image(x2,y1,z2), 'all');
            corner(7) = median(GFP_image(x1,y2,z2), 'all');
            corner(8) = median(GFP_image(x2,y2,z2), 'all');
            min_corner = min(corner);
            
            % Measure the neurons GFP channel.
            ns = neurons.neurons;
            GFP_image = squeeze(data_zscored(:,:,:,GFP_i));
            GFP_colors = nan(length(ns),1);
            for i=1:length(ns)
                cpatch = Methods.Utils.subcube(GFP_image, ...
                    round(ns(i).position), [1,1,0]);
                GFP_colors(i) = median(reshape(cpatch, ...
                    [numel(cpatch)/size(cpatch, 4), size(cpatch, 4)]));
            end
            
            % Measure the neurons to determine an appropriate GFP Otsu threshold.
            colors = neurons.get_colors_readout();
            [otsu_thresh, otsu_score] = graythresh(colors(:,GFP_i));
            otsu_thresh = otsu_thresh * max(colors(:,GFP_i));
            
            % Write the CSV titles & row data to a file
            aligned_xyzRGBs = neurons.get_aligned_xyzRGBs();
            
            % Open the file.
            fileID = fopen(csvfile, 'w');
            
            % Write the background and Otsu thresholds.
            fprintf(fileID, 'Atlas Version,,Background Threshold,,Otsu Threshold,Otsu Score\n');
            fprintf(fileID, '%f,,%f,,%f,%f\n\n', neurons.atlas_version, min_corner, otsu_thresh, otsu_score);
            
            % Determine the header output.
            id_str = 'User ID,User Confidence,Auto ID,Auto Confidence,,';
            real_position_str = 'Real X (um),Real Y (um),Real Z (um),,';
            real_color_str = 'Real Red,Real Green,Real Blue,Real White,,';
            aligned_position_str = [];
            aligned_color_str = [];
            if ~isempty(aligned_xyzRGBs)
                aligned_position_str = 'Aligned X (um),Aligned Y (um),Aligned Z (um),,';
                aligned_color_str = 'Aligned Red,Aligned Green,Aligned Blue,,';
            end
            GFP_str = 'GFP\n';
            out_str = [id_str, real_position_str, real_color_str, ...
                aligned_position_str, aligned_color_str, GFP_str];
            
            % Determine the data output.
            id_fmt = '%s,%f,%s,%f,,';
            real_position_fmt = '%f,%f,%f,,';
            real_color_fmt = '%f,%f,%f,%f,,';
            aligned_position_fmt = [];
            aligned_color_fmt = [];
            if ~isempty(aligned_xyzRGBs)
                aligned_position_fmt = '%f,%f,%f,,';
                aligned_color_fmt = '%f,%f,%f,,';
            end
            GFP_fmt = '%f\n';
            out_fmt = [id_fmt, real_position_fmt, real_color_fmt, ...
                aligned_position_fmt, aligned_color_fmt, GFP_fmt];
                        
            % Write the neurons.
            ns = neurons.neurons;
            fprintf(fileID, out_str);
            for i = 1:length(ns)
                n = ns(i);
                pos = n.position .* um_scale;
                
                % Write the real data only.
                if isempty(aligned_xyzRGBs)
                    fprintf(fileID, out_fmt, ...
                        n.annotation, n.annotation_confidence, ...
                        n.probabilistic_ids{1}, n.probabilistic_probs(1), ...
                        pos(1), pos(2), pos(3), ...
                        n.color_readout(1), n.color_readout(2), n.color_readout(3), n.color_readout(4), ...
                        GFP_colors(i));
                    
                % Write the real & aligned data.
                else
                    aligned_pos = n.aligned_xyzRGB(1:3) .* um_scale;
                    fprintf(fileID, out_fmt, ...
                        n.annotation, n.annotation_confidence, ...
                        n.probabilistic_ids{1}, n.probabilistic_probs(1), ...
                        pos(1), pos(2), pos(3), ...
                        n.color_readout(1), n.color_readout(2), n.color_readout(3), n.color_readout(4), ...
                        aligned_pos(1), aligned_pos(2), aligned_pos(3), ...
                        n.aligned_xyzRGB(4), n.aligned_xyzRGB(5), n.aligned_xyzRGB(6), ...
                        GFP_colors(i));
                end
            end
            
            % Done.
            fclose(fileID);
        end
    end
end
