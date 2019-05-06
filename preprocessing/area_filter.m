function volume = area_filter(volume, point, threshold, dim)
% Heuristic filter for removing objects in a multi-color volume.
%
% Amin Nejat
    point = round(point);
    
    g = volume(:,:,:,dim);
    g(g < threshold) = 0;
    a = bwlabeln(g,26);
    
    mask = logical(0*a);
    mask(a == a(point(1), point(2), point(3))) = 1;
    mask = repmat(mask, 1,1,1,size(volume,4));
    volume(mask) = 0;
    
end

