function [image, metadata] = imreadND2(filename)
%IMREADND2 Read in a Nikon ND2 image.
%
%   [IMAGE, METADATA] = IMREADCZI(FILENAME)
%
%   Input:
%   filename - the filename of the ND2 image
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
xPixelsI = find(contains(keys, 'Global uiWidth'), 1);
yPixelsI = find(contains(keys, 'Global uiHeight'), 1);
numZSlicesI = find(contains(keys, 'Global uiCount'), 1);
xyScaleI = find(contains(keys, 'Global dCalibration'), 1);
zScaleI = find(contains(keys, 'Global dZStep'), 1);

% Extract the image volume information.
image.pixels = [ ...
    values{xPixelsI}; ...
    values{yPixelsI}; ...
    values{numZSlicesI}];
image.scale = [ ...
    values{xyScaleI}; ...
    values{xyScaleI}; ...
    values{zScaleI}];

% Initialize the channel information.
numChannelsI = find(contains(keys, 'Global uiSampleCount'), 1);
channelsI = find(contains(keys, 'Global Name #'));

% Extract the channel information.
numChannels = values{numChannelsI};
channelKeys = keys(channelsI);
channelValues = values(channelsI);
image.channels = cell(numChannels,1);
for i = 1:numChannels
    
    % Get the channel name.
    channelI = find(endsWith(channelKeys, num2str(i + 1)), 1);
    image.channels{i} = channelValues{channelI};
end

% Nikon format is RGBW?
% Note: can't find color info in the metafile, assumming the order is RGBW.
image.colors = nan(numChannels,3);
image.colors(1,:) = [1,0,0];
image.colors(2,:) = [0,1,0];
image.colors(3,:) = [0,0,1];
for i = 4:numChannels
    image.colors(i,:) = [1,1,1];
end
image.dicChannel = nan;

% Initialize the image excitation/emission information.
lasersI = find(contains(keys, ...
    'Global m_uiMultiLaserLineWavelength0-0'));
emissionsI = find(contains(keys, ...
    'Global Lambda 10-B, Filter #'));

% Extract the image excitation/emission information.
laserKeys = keys(lasersI);
laserValues = values(lasersI);
emissionKeys = keys(emissionsI);
emissionValues = values(emissionsI);
image.lasers = nan(numChannels,1);
image.emissions = cell(numChannels,1);
for i=1:numChannels
    
    % Get the laser excitation.
    laserI = find(endsWith(laserKeys, num2str(i-1)), 1);
    image.lasers(i) = laserValues{laserI};
    
    % Get the emission filter.
    emissionI = find(endsWith(emissionKeys, num2str(i)), 1);
    image.emissions{i} = emissionValues{emissionI};  
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
