classdef Illustration
    %ILLUSTRATION Methods to illustrate NeuroPAL images.
    
    % Illustration constants.
    properties (Constant, Access = public)
        font_name = 'Arial';
        font_weight = 'bold';
    end
    
    
   % Public methods.
    methods (Static)
        function saveIDImage(file, data, neurons, um_scale, ...
                start_z, end_z, z_MIP, image_size, font_size)
            %SAVEIDIMAGE save an image with neuron IDs to a PDF.
            
            % Initialize.
            import Output.*;
            
            % Sanitize the input values.
            max_z = size(data,3);
            if start_z < 1 || start_z > max_z
                start_z = 1;
            end
            if end_z < 1 || end_z > max_z
                end_z = max_z;
            end
            if z_MIP < 1 || z_MIP > max_z
                z_MIP = 1;
            end
            if start_z > end_z
                tmp = start_z;
                start_z = end_z;
                end_z = tmp;
            end
            if z_MIP > (end_z - start_z + 1)
                z_MIP = end_z - start_z + 1;
            end
            if image_size <= 0
                image_size = 1;
            end
            if font_size <= 0
                font_size = 10;
            end
            
            % Get the the screen size.
            screen_size = get(0, 'Screensize');
            
            % Compute the scale bar.
            max_y = size(data,1);
            scale_bar_width = 1;
            scale_bar_size = 10;
            scale_bar_str = [num2str(scale_bar_size) ' \mum'];
            scale_bar_off = [10,10];
            scale_bar_len = (scale_bar_size / um_scale(1)) * image_size;
            scale_bar_x = [0, scale_bar_len] + scale_bar_off(1);
            scale_bar_y = max_y * image_size - [scale_bar_off(2),scale_bar_off(2)];
            scale_bar_str_x = sum(scale_bar_x)/2;
            scale_bar_str_y = sum(scale_bar_y)/2;
            
            % Get the neuron info.
            neuron_name = arrayfun(@(x) x.annotation, neurons, 'UniformOutput' , false);
            neuron_pos = vertcat(neurons.position);
            neuron_pos(:,1:2) = neuron_pos(:,1:2) * image_size;
            neuron_conf = arrayfun(@(x) x.annotation_confidence, neurons);
            neuron_on = arrayfun(@(x) x.is_annotation_on, neurons);
            
            % Save the PDFs.
            page = 1;
            for i = start_z:z_MIP:end_z
                
                % Stay in bounds.
                if (i + z_MIP - 1) > end_z
                    z_MIP = end_z - i + 1;
                end
                
                % Draw the image.
                fig = figure('Visible', 'off', 'NumberTitle', 'off', 'Name', file);
                fig.Position(1:2) = [0,screen_size(4)];
                fig.Position(3:4) = [size(data,2),size(data,1)] * image_size;
                z_slice = squeeze(max(data(:,:,i:(i+z_MIP-1),:),[],3));
                if image_size ~= 1
                    z_slice = imresize(z_slice, image_size);
                end
                image_ax = image(z_slice);
                image_ax.Parent.Box = 'off';
                image_ax.Parent.YTick = [];
                image_ax.Parent.XTick = [];
                %image_ax.Parent.PlotBoxAspectRatio = [1,1,1];
                daspect([1,1,1]);
                hold on;
                
                % Find the neurons in these z-slices.
                neuron_z = find(neuron_pos(:,3) >= i & neuron_pos(:,3) <= (i+z_MIP-1));
                
                % Draw the neuron IDs.
                for j = 1:length(neuron_z)
                    k = neuron_z(j);
                    pos = neuron_pos(k,1:2);
                    name = neuron_name{k};
                    
                    % Is the neuron ON/OFF?
                    if ~isnan(neuron_on(k))
                        if neuron_on(k)
                            name = [name '-ON'];
                        else
                            name = [name '-OFF'];
                        end
                    end
                    
                    % Is the user confident about the ID?
                    if neuron_conf(k) < 1
                        name = [name '?'];
                    end
                    
                    %  Draw the ID.
                    text(pos(2), pos(1), name, 'Color', 'w', ...
                        'FontName', Illustration.font_name, ...
                        'FontSize', font_size, ...
                        'FontWeight', Illustration.font_weight, ...
                        'HorizontalAlignment', 'center', ...
                        'VerticalAlignment', 'middle');
                end
                
                % Draw the scale bar.
                plot(scale_bar_x, scale_bar_y, 'w', ...
                    'LineWidth', scale_bar_width);
                text(scale_bar_str_x, scale_bar_str_y, scale_bar_str, ...
                    'Color', 'w', ...
                    'FontName', Illustration.font_name, ...
                    'FontSize', font_size, ...
                    'FontWeight', Illustration.font_weight, ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'bottom');
                
                % Save the file.
                save_file = [file '_' num2str(page) '.pdf'];
                fig.Renderer = 'Painters';
                orient(fig, 'landscape');
                print(fig, '-dpdf', '-fillpage', save_file);
                close(fig);
                page = page + 1;
            end
        end
    end
end
