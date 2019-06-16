function [image, metadata] = imreadCZI(filename)
%IMREADCZI Read in a Zeiss CZI image.
%
%   [IMAGE, METADATA] = IMREADCZI(FILENAME)
%
%   Input:
%   filename - the filename of the CZI image
%
%   Outputs:
%   image - the image, a struct with fields:
%       pixels     = the number of pixels as (x,y,z)
%       scale      = the pixel scale, in meters, as (x,y,z)
%       channels   = the names of the channels
%       colors     = the color for each channel as (R,G,B)
%       dicChannel = the DIC channel number
%       lasers     = the laser wavelength for each channel
%       emissions  = the emssion band for each channel as (min,max)
%       data       = the image data as (x,y,z,channel)
%   metadata - the meta data, a struct with fields:
%       keys      = the meta data keys (the names for the meta data)
%       values    = the meta data values (the values for the meta data)
%       hashtable = a Java Hashtable of keys and their values

% Open the CZI file.
data = bfopen(filename);

% Extract the metadata.
hashtable = data{1,2};
keys = arrayfun(@char, hashtable.keySet.toArray, 'UniformOutput', false);
values = cellfun(@(x) hashtable.get(x), keys, 'UniformOutput', false);

% Organize the metadata.
metadata.keys = keys;
metadata.values = values;
metadata.hashtable = hashtable;

% Initialize the image volume information.
xPixelsI = find(contains(keys, 'Global Information|Image|SizeX #1'), 1);
yPixelsI = find(contains(keys, 'Global Information|Image|SizeY #1'), 1);
xScaleI = find(contains(keys, ...
    'Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1'), 1);
yScaleI = find(contains(keys, ...
    'Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1'), 1);
zScaleI = find(contains(keys, ...
    'Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1'), 1);
numZSlicesI = find(contains(keys, 'Information|Image|SizeZ #1'), 1);

% Extract the image volume information.
image.pixels = [ ...
    str2double(values{xPixelsI}); ...
    str2double(values{yPixelsI}); ...
    str2double(values{numZSlicesI})];
image.scale = [ ...
    str2double(values{xScaleI}); ...
    str2double(values{yScaleI}); ...
    str2double(values{zScaleI})];

% Initialize the channel information.
numChannelsI = find(contains(keys, 'Information|Image|SizeC #1'), 1);
channelsI = find(contains(keys, 'Information|Image|Channel|Name'));
colorsI = find(contains(keys, ...
    'Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Detector|Color'));

% Extract the channel information.
numChannels = floor(str2double(values{numChannelsI}));
channelKeys = keys(channelsI);
channelValues = values(channelsI);
colorKeys = keys(colorsI);
colorValues = values(colorsI);
image.channels = cell(numChannels,1);
image.colors = nan(numChannels,3);
for i=1:numChannels
    
    % Get the channel name.
    channelI = find(endsWith(channelKeys, num2str(i)), 1);
    image.channels{i} = channelValues{channelI};
    
    % Get the channel color.
    colorI = find(endsWith(colorKeys, num2str(i)), 1);
    color = colorValues{colorI};
    image.colors(i,1) = hex2dec(color((end - 5):(end - 4)));
    image.colors(i,2) = hex2dec(color((end - 3):(end - 2)));
    image.colors(i,3) = hex2dec(color((end - 1):end));
end
image.dicChannel = find(contains(image.channels, 'PMT'),1);

% Initialize the image excitation/emission information.
lasersI = find(contains(keys, ...
    'Information|Image|Channel|ExcitationWavelength'));
emissionsI = find(contains(keys, ...
    'Information|Image|Channel|DetectionWavelength|Ranges'));

% Extract the image excitation/emission information.
laserKeys = keys(lasersI);
laserValues = values(lasersI);
emissionKeys = keys(emissionsI);
emissionValues = values(emissionsI);
image.lasers = nan(numChannels,1);
image.emissions = nan(numChannels,2);
for i=1:numChannels
    
    % The DIC channel has no laser nor an emission band.
    if i == image.dicChannel
        continue;
    end
    
    % Skip over the DIC channel.
    bandI = i;
    if i > image.dicChannel
        bandI = i - 1;
    end
    
    % Get the laser excitation.
    laserI = find(endsWith(laserKeys, num2str(bandI)), 1);
    image.lasers(i) = str2double(laserValues{laserI});
    
    % Get the emission band.
    emissionI = find(endsWith(emissionKeys, num2str(bandI)), 1);
    image.emissions(i,:) = sscanf(emissionValues{emissionI}, '%f-%f');  
end

% Organize the image volume.
imageData = data{1,1};
image.data = uint16(nan([image.pixels; numChannels]'));
numC = numChannels;
for i=1:size(imageData,1)
    
    % Assemble the image.
    z = floor((i - 1) / numC) + 1;
    c = mod(i - 1, numC) + 1;
    image.data(:,:,z,c) = imageData{i,1}';
    
    % Debug the image assembly.
    % disp(imageData{i,2});
    % printf('z=%d c=%d', z, c);
end
end
