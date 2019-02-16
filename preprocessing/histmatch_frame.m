function hmvideo = histmatch_frame(video)
% matching the histogram of the intensities for different colors to
% normalize different color channels for better visualization
% Amin Nejat

hmvideo = zeros(size(video));
exprand = exprnd(1, 1, numel(find(video > 0)));

for ch = 1: size(video, 4)
    vol = squeeze(video(:,:,:,ch,:));
    posvol = vol(vol > 0);

    result = vol;
    result(vol > 0) = imhistmatchn(posvol(:)/max(posvol(:)), exprand/max(exprand(:)), 100000);
    hmvideo(:,:,:,ch,:) = reshape(result, size(vol));
end
    
        
end