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