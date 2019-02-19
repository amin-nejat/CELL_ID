function thresh = threshold_frame(video, method, sensitivity)
% Thresholding multi-color volume
%
% Amin Nejat

thresh = 0*video;

for ch = 1: size(video, 4)
    vol = squeeze(video(:,:,:,ch,:));
    if strcmp(method, 'adaptive')
        BW = imbinarize(vol/max(vol(:)), 'adaptive', 'Sensitivity', sensitivity);
        vol(~BW) = 0;
        vol = imclearborder(vol);
        thresh(:,:,:,ch,:) = vol;
    elseif strcmp(method, 'global')
        BW = imbinarize(vol/max(vol(:)), 'global');
        vol(~BW) = 0;
        vol = imclearborder(vol);
        thresh(:,:,:,ch,:) = vol;
    elseif strcmp(method, 'otsu')
        vol(vol < graythresh(vol(:))) = 0;
        thresh(:,:,:,ch,:) = vol;
    elseif strcmp(method, 'fixed')
        vol(vol < sensitivity) = 0;
        thresh(:,:,:,ch,:) = vol;
    end
end

end