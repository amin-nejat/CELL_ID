classdef Illustration
    %ILLUSTRATION Methods to illustrate NeuroPAL images.
    
    % Illustration constants.
    properties (Constant, Access = public)
        font_name = 'Arial';
        font_weight = 'bold';
    end
    
    
   % Public methods.
    methods (Static)
        function saveIDImage(file, data, neurons, um_scale)
            %SAVEIDIMAGE save an image with neuron IDs to a PDF.
            
            % Initialize.
            import Output.*;
            
            % What is the z limit?
            max_z = size(data,3);
            max_z_str = num2str(max_z);
            range_z_str = ['(1-' max_z_str ')'];
            range_z_ID_str = ['(0-' max_z_str ')'];
            
            % Worm neurons are ~2-4um in diameter so default to 8um MIP z-slices.
            z_MIP_def = num2str(round(8 / um_scale(3)));
            
            % Prompt the user for the parameters.
            start_z_str = ['Start Z-slice ' range_z_str];
            end_z_str = ['End Z-slice ' range_z_str];
            z_MIP_str = ['Z-slice thickness ' range_z_str];
            z_ID_str = ['Z-slice ID thickness ' range_z_ID_str];
            size_str = 'Image size (default = 1x)';
            font_str = 'Font size';
            prompt = {start_z_str, end_z_str, z_MIP_str, z_ID_str, ...
                size_str, font_str};
            title = 'Save ID Image';
            dims = [1 35];
            definput = {'1', max_z_str, z_MIP_def, '1', '2', '4'};
            answer = inputdlg(prompt, title, dims, definput);
            if isempty(answer)
                return;
            end
            
            % Format the user input.
            start_z = round(str2double(answer{1}));
            end_z = round(str2double(answer{2}));
            z_MIP = round(str2double(answer{3}));
            z_ID_offset = round(str2double(answer{4}));
            image_size = str2double(answer{5});
            font_size = round(str2double(answer{6}));
            
            % Sanitize the input values.
            if start_z < 1 || start_z > max_z
                start_z = 1;
            end
            if end_z < 1 || end_z > max_z
                end_z = max_z;
            end
            if z_MIP < 1 || z_MIP > max_z
                z_MIP = 1;
            end
            if z_ID_offset < 0 || z_ID_offset > max_z
                z_ID_offset = 1;
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
            
            % Take a minimal patch around the neuron center.
            % Note: we need to walk a thin line of being robust against
            % off-center dots, while not violating neighboring neurons.
            cube_size = [1,1,1];
                
            % Get the neuron info.
            neuron_pos = vertcat(neurons.position);
            if ~isempty(neuron_pos)
                
                % Compute the neuron brightness for overlaying text.
                neuron_brightness_thresh = 0.7;
                image_brightness = data(:,:,:,2); % use the green channel
                neuron_brightness = arrayfun(@(x) median( ...
                    Methods.Utils.subcube(image_brightness, ...
                    round(neuron_pos(x,:)), cube_size), 'all'), ...
                    1:size(neuron_pos,1));
                
                % Determine the neuron names, image coordinates, ID
                % confidence, & ON/OFF annotations.
                neuron_name = arrayfun(@(x) x.annotation, neurons, 'UniformOutput' , false);
                neuron_pos(:,1:2) = neuron_pos(:,1:2) * image_size;
                neuron_conf = arrayfun(@(x) x.annotation_confidence, neurons);
                neuron_on = arrayfun(@(x) x.is_annotation_on, neurons);
            end
            
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
                neuron_z = [];
                if z_ID_offset > 0 && ~isempty(neuron_pos)
                    min_z_ID = i - z_ID_offset + 1;
                    max_z_ID = i + z_MIP + z_ID_offset - 2;
                    neuron_z = find(neuron_pos(:,3) >= min_z_ID & ...
                        neuron_pos(:,3) <= max_z_ID);
                end
                
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
                    
                    % Are we writing on a bright green background?
                    text_color = 'w';
                    if neuron_brightness(k) > neuron_brightness_thresh
                        text_color = 'k';
                    end
                    
                    % Draw the ID.
                    text(pos(2), pos(1), name, 'Color', text_color, ...
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
