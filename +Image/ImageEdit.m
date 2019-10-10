classdef ImageEdit
    %IMAGEEDIT Image editing methods.
    
    methods (Static)
        function image = convert2uint16(image)
            %CONVERT2UINT16 Convert the image to unit16.
            image = image - min(image(:));
            image = uint16(double(intmax('uint16')) * double(image) ...
                / double(max(image(:))));
        end
        
        function image = convert2double01(image)
            %CONVERT2DOUBLE01 Convert the image to [0,1].
            image = double(image);
            image = image - min(image(:));
            image = image / max(image(:));
        end
        
        function zimage = zscore(image)
            %ZSCORE Z-score the image colors.
            
            % Z-score the image.
            zimage = double(image);
            for c = 1:size(image,4)
                data = zimage(:,:,:,c,:);
                zimage(:,:,:,c,:) = (data - mean(data(:))) / std(data(:));
            end
            
            % Convert the image to a viewable format.
            zimage = Image.ImageEdit.convert2double01(zimage);
        end
        
        function bimage = deBleed(bimage, dsts, srcs, scales)
            %DEBLEED De-bleed the destination channel(s) of their source
            % channel(s) bleedthrough using the specified scale(s).
            
            % De-bleed every channel using the same scale.
            if length(scales) == 1
                bimage(:,:,:,dsts) = bimage(:,:,:,dsts) - ...
                    bimage(:,:,:,srcs) * scales;
                
                % De-bleed each channel using its individual scale.
            else
                for i=1:length(scales)
                    bimage(:,:,:,dsts(i)) = bimage(:,:,:,dsts(i)) - ...
                        bimage(:,:,:,srcs(i)) * scales(i);
                end
            end
        end
        
        function aimage = adjust(aimage, lows, highs, gammas)
            %ADJUST Adjust the image histogram and gamma.
            
            % Setup the defaults.
            if isempty(lows)
                lows = 0;
            end
            if isempty(highs)
                highs = 1;
            end
            if isempty(gammas)
                gammas = 1;
            end
            
            % Are we adjusting each color channel separately?
            lows((length(lows)+1):size(aimage,4)) = lows(1);
            highs((length(highs)+1):size(aimage,4)) = highs(1);
            gammas((length(gammas)+1):size(aimage,4)) = gammas(1);
            
            % Adust the color channels.
            for c = 1:size(aimage,4)
                aimage(:,:,:,c) = imadjustn(squeeze(aimage(:,:,:,c)), ...
                    [lows(c) highs(c)], [], gammas(c));
            end
        end
        
        function himage = adjustHistogram(himage, lows, highs)
            %ADJUSTHISTOGRAM Adjust the image histogram.
            himage = Image.ImageEdit.adjust(himage, lows, highs, 1);
        end
        
        function gimage = adjustGamma(gimage, gammas)
            %ADJUSTGAMMA Adjust the image gamma.
            gimage = Image.ImageEdit.adjust(gimage, [], [], gammas);
        end

        function himage = deHaze(himage, RGB, varargin)
            %DEHAZE Remove haze.
            
            % Saturate the image.
            saturation = 1.5;
            if ~isempty(varargin)
                saturation = varargin{1};
            end
            
            % De-haze the image.
            for z = 1:size(himage, 3)
                
                % Convert to CIELAB color space for lightness manipulation.
                h = rgb2lab(squeeze(himage(:,:,z,RGB)));
                
                % Reduce haze in the lightness channel.
                lum = h(:,:,1) ./ 100;
                lum = imreducehaze(lum, 'ContrastEnhancement', 'none');
                h2 = lum;
                h2(:,:,1) = h2 .* 100;
                
                % Saturate the color channels.
                h2(:,:,2:3) = h(:,:,2:3) * saturation;
                
                % Convert back to RGB color space.
                h2 = lab2rgb(h2);
                himage(:,:,z,RGB) = h2;
            end
            
            % Convert the image to a viewable format.
            himage = Image.ImageEdit.convert2uint16(himage);
        end
        
        function mimage = morphOpen(mimage, varargin)
            %MORPHOPEN Morphologically open the image, removing small
            % background objects (noise).
            
            % How big are the objects we want to remove?
            strel_size = 1;
            if ~isempty(varargin)
                strel_size = varargin{1};
            end
            
            % Create the structuring element.
            % Note: most microscope images are anisotropic wherein Z is
            % lower resolution than X and Y. Therefore, a sphere is an
            % inappropriate structuring element.
            strel_type = 'disk'; % 'disk'
            se = strel(strel_type, strel_size);
            
            % Morphologically open the image.
            for c = 1:size(mimage, 4)
                I = squeeze(mimage(:,:,:,c));
                I = imopen(I, se);
                mimage(:,:,:,c) = I;
            end
        end
        
        function nimage = deNoise(nimage)
            %DENOISE Remove noise.
            for c = 1:size(nimage, 4)
                for z = 1:size(nimage, 3)
                    I = squeeze(nimage(:,:,z,c));
                    nimage(:,:,z,c) = imdiffusefilt(I);
                end
            end
        end
    end
end
