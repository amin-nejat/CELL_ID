function [video, shapes] = sp2video(sp, sz, nsz, trunc)
% Reconstruction of the MoG, this function converts a superpixel data
% structure to a video.
%
% Amin Nejat

video = zeros(sz);

shapes = cell(sz(5), size(sp(1).mean,1)); 

for t = 1: sz(5)
    for neuron = 1: size(sp(1).mean, 1)
        
        if size(sp(t).baseline, 1) > 1
            baseline_mean = sp(t).baseline(neuron, :);
        else
            baseline_mean = sp(t).baseline;
        end
        
        baseline_std = zeros(length(baseline_mean));

        abs_loc = sp(t).mean(neuron, :);
        cov = squeeze(sp(t).cov(neuron, :, :));
    
        color = sp(t).color(neuron, :);


    
        
        simul = simulate_gaussian(2*nsz-1, ... % size
            nsz, ... % center
            cov, ... % cov
            color, ... % colors
            0*baseline_mean, ... % baseline mean
            0*baseline_std, ... % baseline std
            trunc); % truncation
        shapes{t, neuron} = simul;
        video(:, :, :, :, t) = video(:, :, :, :, t) + placement(sz(1:3), round(abs_loc(1:3)), simul);
    end
    
end

end