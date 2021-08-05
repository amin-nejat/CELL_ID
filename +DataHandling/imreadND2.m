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

% Open the ND2 file.
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
xPixelsI = find(contains(keys, 'Global uiWidth') ...
    & ~contains(keys, 'Bytes'), 1);
yPixelsI = find(contains(keys, 'Global uiHeight') ...
    & ~contains(keys, 'Bytes'), 1);
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
numChannelsI = find(contains(keys, 'Global Number of Picture Planes'), 1);
channelsI = find(contains(keys, 'Global Name #'));

% Extract the channel information.
numChannels = round(str2double(values{numChannelsI}));
channelKeys = keys(channelsI);
channelValues = values(channelsI);
image.channels = cell(numChannels,1);
for i = 1:numChannels
    
    % Get the channel name.
    channelI = find(endsWith(channelKeys, num2str(i + 1)), 1);
    image.channels{i} = channelValues{channelI};
end

% Default to BGRW.
% Note: can't find color info in the metafile.
red = [1,0,0];
green = [0,1,0];
blue = [0,0,1];
white = [1,1,1];
image.colors = nan(numChannels,3);
image.colors(1,:) = blue;
image.colors(2,:) = green;
image.colors(3,:) = red;
for i = 4:numChannels
    image.colors(i,:) = white;
end
image.dicChannel = nan;

% Try using the channel names to determine their colors.
for i = 1:length(image.channels)
    wavelength = sscanf(image.channels{i}, '%f');
    
    % DIC.
    if isempty(wavelength) || isnan(wavelength) || wavelength < 350
        image.dicChannel = i;
        image.colors(i,:) = white;
    elseif wavelength < 440
        image.colors(i,:) = blue;
    elseif wavelength < 530
        image.colors(i,:) = green;
    else
        image.colors(i,:) = red;
    end
end

% Initialize the image excitation/emission information.
% lasersI = find(contains(keys, ...
%     'Global m_uiMultiLaserLineWavelength0-0'));
% emissionsI = find(contains(keys, 'Lambda') & contains(keys, 'Filter'));

% Extract the image excitation/emission information.
% laserKeys = keys(lasersI);
% laserValues = values(lasersI);
% emissionKeys = keys(emissionsI);
% emissionValues = values(emissionsI);
% image.lasers = nan(numChannels,1);
% image.emissions = cell(numChannels,1);
% for i=1:numChannels
%     
%     % Get the laser excitation.
%     laserI = find(endsWith(laserKeys, num2str(i-1)), 1);
%     image.lasers(i) = laserValues{laserI};
%     
%     % Get the emission filter.
%     emissionI = find(endsWith(emissionKeys, num2str(i)), 1);
%     image.emissions{i} = emissionValues{emissionI};  
% end

% Use the lasers to determine the color channels.
% Note: assume initial color channel assignments take precedence.
% for i=flip(1:size(image.lasers,1))
%     % mTagBFP2
%     if image.lasers(i) < 440
%         image.colors(i,:) = blue;
%     % CyOFP (or GFP).
%     elseif image.lasers(i) < 530
%         image.colors(i,:) = green;
%     % TagRFP-T
%     elseif image.lasers(i) < 570
%         image.colors(i,:) = white;
%     % mNeptune2.5
%     elseif image.lasers(i) < 650
%         image.colors(i,:) = red;
%     % Anything else.
%     else
%         image.colors(i,:) = white;
%     end
% end

% Organize the image volume.
%numC = numChannels;
imageData = data{1,1};
image.data = uint16(nan([image.pixels; numChannels]'));
for i=1:size(imageData,1)
    
    % Get the image plane data.
    dataStrs = split(imageData{i,2}, ';');
    zStr = strtrim(dataStrs{end-1});
    cStr = strtrim(dataStrs{end});
    
    % Assemble the image.
    %z = floor((i - 1) / numC) + 1;
    %c = mod(i - 1, numC) + 1;
    %image.data(:,:,z,c) = imageData{i,1}';
    z = sscanf(zStr,'Z=%f');
    c = sscanf(cStr,'C=%f');
    if isempty(z) || isnan(z)
        z = 1;
    end
    if isempty(c) || isnan(c)
        c = 1;
    end
    image.data(:,:,z,c) = imageData{i,1}';
    
    % Debug the image assembly.
    % disp(imageData{i,2});
    % printf('z=%d c=%d', z, c);
end
end
