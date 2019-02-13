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