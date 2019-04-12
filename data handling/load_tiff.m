function [bbVol, info] = load_tiff(imagefile)
% Load tiff files in MATLAB. Requires mij.jar and Fiji software.
%
% Amin Nejat
params = get_params();
addpath(genpath(params.fiji_lib_folder));
Miji(false);
MIJ.run('Open...', ['path=[', imagefile, ']']);
thisVol = MIJ.getCurrentImage;


[~, ~, ext] = fileparts(imagefile);
ext = ext(1: 4);
switch ext
    case '.czi'
        iminfo = czifinfo(imagefile);
        desc = iminfo.metadataXML;

        names = regexp(desc, '<DimensionZ>(?<slices>\d+)</DimensionZ>', 'names');
        try slices = str2num(names.slices); catch slices = 0; end

        names = regexp(desc, '<DimensionT>(?<frames>\d+)</DimensionT>', 'names');
        try frames = str2num(names.frames); catch frames = 1; end
        frames = 1;

        names = regexp(desc, '<DimensionX>(?<x>\d+)</DimensionX>', 'names');
        try x = str2num(names.x); catch x = 0; end

        names = regexp(desc, '<DimensionY>(?<y>\d+)</DimensionY>', 'names');
        try y = str2num(names.y); catch y = 0; end

        ch_count = size(thisVol, 3)/(slices*frames);

        info = iminfo;

    case '.tif'
        info    = imfinfo(imagefile);

        try desc = info(1).ImageDescription;
            names = regexp(desc, 'SizeC="(?<channels>\d+)"\s*', 'names');
            try ch_count = str2num(names.channels); catch ch_count = info(1).SamplesPerPixel; end

            names = regexp(desc, 'slices=(?<slices>\d+)\s*', 'names');
            try slices = str2num(names.slices); catch slices = 0; end

            names = regexp(desc, 'frames=(?<frames>\d+\s*)', 'names');
            try frames = str2num(names.frames); catch frames = 1; end
        catch
            ch_count = 1; slices = size(thisVol, 3); frames = 1; 
        end

        info = info(1);
end



if slices == 0
    slices = numel(thisVol)/prod([size(thisVol, 1), size(thisVol, 2), ch_count, frames]);
end

bbVol = reshape(thisVol, [size(thisVol, 1), size(thisVol, 2), slices, ch_count, frames]);
% bbVol = permute(bbVol, [1, 2, 4, 3, 5]);

end