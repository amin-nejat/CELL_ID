classdef NeuroPALImage
    %NEUROPALIMAGE Convert various image formats to a NeuroPAL format.
    %
    %   NeuroPAL files contain 4 variable:
    %      data = the image data (x,y,z,c)
    %      info = the image information
    %           scale = pixel scale in microns (x,y,z)
    %           RGBW = the (R,G,B,W) color channel indices (nan = no data)
    %           GFP = the GFP color channel index(s) (can be empty)
    %           DIC = the DIC channel index (can be empty)
    %           gamma = the gamma correction for the image
    %      prefs = the user preferences
    %           RGBW = the (R,G,B,W) color channel indices (nan = no data)
    %           GFP = the GFP color channel index(s) (can be empty)
    %           DIC = the DIC channel index (can be empty)
    %           gamma = the gamma correction for the image
    %           rotate.horizontal = rotate horizontal?
    %           rotate.vertical = rotate vertical?
    %           body_part = 'Head', 'Midbody', 'Tail'
    %      neurons = the neurons in the image

    
    %% Public methods.
    methods (Static)
        function [data, info, prefs, neurons, np_file] = open(filename)
            %OPEN Open an image in NeuroPAL format.
            %
            %   data = the image data
            %   info = the image information
            %   prefs = the user preferences
            %   neurons = the neurons in the image
            %   filename = the NeuroPAL format filename
            
            % Initialize the return values.
            %data = [];
            %info = [];
            %prefs = [];
            %neurons = [];
            
            % Is the user accidentally trying to open the ID file?
            id_file_ext = '_ID.mat';
            if endsWith(filename, id_file_ext)
                filename = strrep(filename, id_file_ext, '.mat');
            end
            
            % Get the file extension.
            [~, ~, ext] = fileparts(filename);
            if isempty(ext)
                error('Unknown image format: "%s"', filename);
            end
            ext = lower(ext);
            
            % Determine the NeuroPAL filename.
            np_file = strrep(filename, ext, '.mat');
            
            % Is the file already in NeuroPAL format?
            if ~exist(np_file,'file')
                switch lower(ext)
                    case '.mat' % NeuroPAL format
                        error('File not found: "%s"', filename);
                    case '.czi' % Zeiss format
                        NeuroPALImage.convertCZI(filename);
                    case '.nd2' % Nikon format
                        % TODO
                    case {'.tif','.tiff'} % TIFF format
                        % TODO
                    otherwise % Unknown format
                        error('Unknown image format: "%s"', filename);
                end
            end
            
            % Open the file.
            np_data = load(np_file);
            if ~isfield(np_data, 'data') || ...
                    ~isfield(np_data, 'info') || ...
                    ~isfield(np_data, 'prefs') || ...
                    ~isfield(np_data, 'neurons')
                error('Misformatted NeuroPAL file: "%s"', filename);
            end
            
            % Setup the file contents.
            data = np_data.data;
            info = np_data.info;
            prefs = np_data.prefs;
            neurons = np_data.neurons;
        end
    end
    
    
    %% Private methods.
    methods (Static, Access = private)
        function np_file = convertCZI(czi_file)
            %CONVERTCZI Convert a CZI file to NeuroPAL format.
            %
            % czi_file = the CZI file to convert
            % np_file = the NeuroPAL format file
            
            % Open the file.
            np_file = [];
            try
                [image_data, meta_data] = imreadCZI(czi_file);
            catch
                warning('Cannot read: "%s"', czi_file);
                return;
            end
            
            % Check the image orientation.
            data = image_data.data;
            data_order = 1:ndims(data);
            if size(data,1) > size(data,2)
                
                % Fix the orientation.
                data_order(1) = 2;
                data_order(2) = 1;
                data = permute(data, data_order);
                
                % Reorder the image scale.
                scale = image_data.scale;
                image_data.scale(1) = scale(2);
                image_data.scale(2) = scale(1);
            end
            
            % Setup the NP file data.
            info.scale = image_data.scale * 1000000; % convert to microns
            info.DIC = image_data.dicChannel;
            
            % Determine the color channels.
            colors = image_data.colors;
            colors = round(colors/max(colors(:)));
            info.RGBW = nan(4,1);
            info.GFP = nan;
            for i = 1:size(colors,1)
                switch char(colors(i,:))
                    case [1,0,0] % red
                        info.RGBW(1) = i;
                    case [0,1,0] % green
                        info.RGBW(2) = i;
                    case [0,0,1] % blue
                        info.RGBW(3) = i;
                    case [1,1,1] % white
                        if i ~= info.DIC
                            info.RGBW(4) = i;
                        end
                    otherwise % GFP
                        info.GFP = i;
                end
            end
            
            % Did we find the GFP channel?
            if isnan(info.GFP) && size(colors,1) > 4
                
                % Assume the first unused channel is GFP.
                unused = setdiff(1:size(colors,1), info.RGBW);
                info.GFP = unused(1);
            end
            
            % Determine the gamma.
            info.gamma = 1;
            keys = lower(meta_data.keys);
            gamma_i = find(contains(keys, 'gamma'),1);
            if ~isempty(gamma_i)
                info.gamma = str2double(meta_data.values(gamma_i));
            end
            
            % Initialize the user preferences.
            prefs.RGBW = info.RGBW;
            prefs.DIC = info.DIC;
            prefs.GFP = info.GFP;
            prefs.gamma = info.gamma;
            prefs.rotate.horizontal = false;
            prefs.rotate.vertical = false;
            prefs.body_part = [];
            
            % Save the CZI file to our MAT file format.
            np_file = strrep(czi_file, 'czi', 'mat');
            neurons = [];
            save(np_file, 'data', 'info', 'prefs', 'neurons');
        end
    end
end

