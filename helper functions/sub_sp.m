function subsp = sub_sp(sp, subset)
% Calculating a subset of superpixels.
%
% Amin Nejat

    subsp = sp;
        
    for t = 1: length(sp)
        subsp(t).mean = sp(t).mean(subset, :);
        subsp(t).color = sp(t).color(subset, :);
        subsp(t).cov = sp(t).cov(subset, :, :);
        try
            subsp(t).baseline = sp(t).baseline(subset, :);
        catch
            subsp(t).baseline = sp(t).baseline;
        end
    end
    
end