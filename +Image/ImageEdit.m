classdef ImageEdit
    %IMAGEEDIT Image editing methods.
    
    methods (Static)
        function zimage = zscore(image)
        %ZSCORE Z-score the image colors.
        zimage = double(image);
        for c = 1:size(image,4)
            data = zimage(:,:,:,c,:);
            zimage(:,:,:,c,:) = (data - mean(data(:))) / std(data(:));
        end
        end
        
        function image = adjust(image, low, high, gamma)
        %ADJUST Adjust the image histogram and gamma.
        hist_thresh = [low high];
        for c = 1:size(image,4)
            image(:,:,:,c) = imadjustn(squeeze(image(:,:,:,c)), ...
                hist_thresh, [], gamma);
        end
        end
        
        function image = adjustHistogram(image, low, high)
        %ADJUSTHISTOGRAM Adjust the image histogram.
            image = adjust(image, low, high, []);
        end
        
        function image = adjustGamma(image, gamma)
        %ADJUSTGAMMA Adjust the image gamma.
            image = adjust(image, [], [], gamma);
        end
    end
end

