classdef ImageEdit
    %IMAGEEDIT Image editing methods.
    
    methods (Static)
        function zimage = zscore(image)
        %ZSCORE Z-score the image colors.
        
        % Z-score the image.
        zimage = double(image);
        for c = 1:size(image,4)
            data = zimage(:,:,:,c,:);
            zimage(:,:,:,c,:) = (data - mean(data(:))) / std(data(:));
        end
        
        % Convert the image to a viewable format.
        zimage = uint16(double(intmax('uint16')) * double(zimage) ...
            / double(max(zimage(:))));
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
        
        function bimage = deBleed(bimage, dsts, srcs, scales)
            %DEBLEED De-bleed the destination channel(s) of their source
            % channel(s) bleedthrough using the specified scale(s).
            bimage(:,:,:,dsts) = bimage(:,:,:,dsts) - ...
                bimage(:,:,:,srcs) .* scales;
        end
    end
end

