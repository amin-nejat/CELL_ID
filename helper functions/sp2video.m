function [video, shapes] = sp2video(sp, sz, nsz, trunc)

shapes = cell(sz(5), size(sp(1).mean,1));
video = zeros(sz);
for t = 1: sz(5)
    for n=1:size(sp(t).mean,1)
        n
        simul = simulate_gaussian(2*nsz+1, ... % size
                nsz+1+sp(t).mean(n,:)-round(sp(t).mean(n,:)), ... % center
                squeeze(sp(t).cov(n,:,:)), ... % cov
                sp(t).color(n,:), ... % colors
                zeros(size(sp(t).color(n,:))), ... % baseline mean
                zeros(size(sp(t).color(n,:))), ... % baseline std
                trunc); % truncation
        shapes{t,n} = simul;
        video(:,:,:,:,t) = superpose(video(:,:,:,:,t), round(sp(t).mean(n,:)), simul);
        
    end
end

end