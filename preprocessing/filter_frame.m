function rho = filter_frame(frame, filter)
% applying a filter to multi-color RGB image
% Amin Nejat

rho = 0*frame;
for ch = 1: size(frame, 4)
    D = squeeze(frame(:, :, :, ch));
    rho(:, :, :, ch) =  imfilter(D, filter, 'same');
end

end

