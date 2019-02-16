function areaopen = rm_smallobj_frame(volume, minsize)
% Heuristic filter for removing small objects in a multi-color volume.
%
% Amin Nejat

    areaopen = 0*volume;
    for ch = 1: size(areaopen, 4)
        tmp = squeeze(volume(:, :, :, ch));
        bwao = bwareaopen(tmp, minsize, 26);
        tmp(~bwao) = 0;
        areaopen(:, :, :, ch) = tmp;
    end
end

