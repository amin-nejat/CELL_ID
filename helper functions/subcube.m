function patch = subcube(cube, loc, center)
% Grabbing a patch of the image in a given location.
%
% Amin Nejat

    sz = [size(cube, 1), size(cube, 2), size(cube, 3)];
    
    rel = round(center);
    reu = round(center);
    
    rel(loc - center - 1 < 0) = loc(loc - center - 1 < 0) - 1;
    reu(loc + center - sz > 0) = sz(loc + center - sz > 0) - loc(loc + center - sz > 0);
                
    patch = cube(loc(1)-rel(1): loc(1)+reu(1), ...
                 loc(2)-rel(2): loc(2)+reu(2), ...
                 loc(3)-rel(3): loc(3)+reu(3), :);
    newcenter = [size(patch, 1), size(patch, 2), size(patch, 3)];
    
    if any(newcenter(1: 3) ~= 2*round(center)+1)
        pre = round(center) - rel;
        post = round(center) - reu;
        
        patch = padarray(padarray(patch, pre, 'pre'), post, 'post');
    end
end