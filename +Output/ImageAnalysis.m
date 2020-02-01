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

            % Get the GFP channel.
            % Note: we translate the z-score to >= 0. Negative GFP
            % intensities confuse users.
            GFP_i = prefs.GFP;
            GFP_image = squeeze(data_zscored(:,:,:,GFP_i));
            GFP_image = GFP_image - min(GFP_image(:));
            
            % Measure the 8 image corners to determine an appropriate
            % GFP background threshold.
            GFP_bg_size = Output.ImageAnalysis.GFP_bg_size;
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
            
            % Take a minimal patch around the neuron center.
            % Note: we need to walk a thin line of being robust against
            % off-center dots, while not violating neighboring neurons.
            cube_size = [1,1,1];
            % Take a 1 micron radius around the neuron centers.
            % Note: users are way too sloppy for us to use this :(
            %cube_size = round([1,1,1] ./ um_scale');
            
            % Measure the GFP channel for the neurons.
            intensity_prctile = 75;
            ns = neurons.neurons;
            GFP_colors = nan(length(ns),1);
            for i=1:length(ns)
                cpatch = Methods.Utils.subcube(GFP_image, ...
                    round(ns(i).position), cube_size);
                GFP_colors(i) = mean(cpatch(cpatch > prctile(cpatch(:), ...
                    intensity_prctile)));
            end
            
            % Compute an appropriate GFP Otsu threshold.
            [otsu_thresh, otsu_score] = graythresh(GFP_colors);
            otsu_thresh = otsu_thresh * max(GFP_colors);
            
            % Compute an appropriate GFP linear change threshold.
            GFP_sorted = sort(GFP_colors);
            [change_point, change_residual] = findchangepts(GFP_sorted, ...
                'MaxNumChanges', 1, 'Statistic', 'linear');
            change_thresh = nan;
            change_i = round(change_point);
            if change_i > 1 && change_i < length(GFP_sorted)
                change_thresh = GFP_sorted(change_i);
            end

            % Get the aligned neuron data.
            aligned_xyzRGBs = neurons.get_aligned_xyzRGBs();
            %if ~isempty(aligned_xyzRGBs)
                
                % Get the real neuron colors.
                %neuron_RGBWs = neurons.get_colors_readout();
                %neuron_RGBs = neuron_RGBWs(:,[1 2 3]);
            
                % Find transformation between original and aligned data.
                %beta_col = linsolve([neuron_RGBs ones(size(neuron_RGBs,1),1)],...
                %    [aligned_xyzRGBs(:,4:end) ones(size(neuron_RGBs,1),1)]);
                %aligned_colors_rgb = [neuron_RGBWs(:,[1,2,3]) ones(size(neuron_RGBWs,1),1)]*beta_col;
                %aligned_colors_rab = [neuron_RGBWs(:,[1,4,3]) ones(size(neuron_RGBWs,1),1)]*beta_col;
                %aligned_colors = [aligned_colors_rgb(:,1:3), aligned_colors_rab(:,2)];
            %end
            
            % Open the file.
            fileID = fopen(csvfile, 'w');
            
            % Write the background and Otsu thresholds.
            fprintf(fileID, ['Atlas Version,,Background Threshold,,' ...
                'Otsu Threshold,Otsu Score,,' ...
                'Linear Change Point,Change Residual\n']);
            fprintf(fileID, '%f,,%f,,%f,%f,,%f,%f\n\n', neurons.atlas_version, ...
                min_corner, otsu_thresh, otsu_score, ...
                change_thresh, change_residual);
            
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
            
            % Sort the neurons by position.
            % Note: x & y are reversed.
            positions = neurons.get_positions();
            positions = positions(:,[2,1,3]);
            [~, sort_i] = sortrows(positions);
            ns = ns(sort_i);
            GFP_colors = GFP_colors(sort_i);
            
            % Write the neurons.
            um_scale = um_scale';
            fprintf(fileID, out_str);
            for i = 1:length(ns)
                n = ns(i);
                pos = n.position .* um_scale;
                
                % Determine the auto IDs.
                probabilistic_id = [];
                probabilistic_prob = [];
                if ~isempty(n.probabilistic_ids)
                    probabilistic_id = n.probabilistic_ids{1};
                    probabilistic_prob = n.probabilistic_probs(1);
                end
                
                % Write the real data only.
                if isempty(aligned_xyzRGBs)
                    fprintf(fileID, out_fmt, ...
                        n.annotation, n.annotation_confidence, ...
                        probabilistic_id, probabilistic_prob, ...
                        pos(2), pos(1), pos(3), ...
                        n.color_readout(1), n.color_readout(2), n.color_readout(3), n.color_readout(4), ...
                        GFP_colors(i));
                    
                % Write the real & aligned data.
                else
                    
                    % Does this neuron have aligned data?
                    aligned_pos = nan(1,3);
                    aligned_RGB = nan(1,3);
                    if ~isempty(n.aligned_xyzRGB)
                        aligned_pos = n.aligned_xyzRGB(1:3) .* um_scale;
                        aligned_RGB = n.aligned_xyzRGB(4:6);
                    end
                    
                    % Write the real & aligned data.
                    fprintf(fileID, out_fmt, ...
                        n.annotation, n.annotation_confidence, ...
                        probabilistic_id, probabilistic_prob, ...
                        pos(2), pos(1), pos(3), ...
                        n.color_readout(1), n.color_readout(2), ...,
                        n.color_readout(3), n.color_readout(4), ...
                        aligned_pos(2), aligned_pos(1), aligned_pos(3), ...
                        aligned_RGB(1), aligned_RGB(2), aligned_RGB(3), ...
                        GFP_colors(i));
                end
            end
            
            % Done.
            fclose(fileID);
        end
    end
end
