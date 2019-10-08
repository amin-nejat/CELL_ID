function [image, metadata] = imreadLif(filename)
%IMREADANY Read .lif format (all purpose reader).
%
%   [IMAGE, METADATA] = IMREADANY(FILENAME)
%
%   Input:
%   filename - the filename of the lif image
%   idx - index of the series (0 to n-1)
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



% Initialize the image data.

image = [];
metadata = [];

% Open the file.

% create reader
reader=bfGetReader(filename); %create a reader using the .lif file name as input
glob=reader.getGlobalMetadata();
ser=reader.getSeriesMetadata();
javaMethod('merge', 'loci.formats.MetadataTools', ...
    glob, ser, 'Global ');


omeMeta=reader.getMetadataStore();
nSeries=reader.getSeriesCount();

% Get dimensions of the dataset
sx=reader.getSizeX();
sy=reader.getSizeY();
sz=reader.getSizeZ();
sC=reader.getSizeC();



DatInf={'1','405','488','561','630','0'};

prompt=cell(1,sC+1);
prompt{1}=['Series to load (out of ' int2str(nSeries) ')'];
for i=1:sC
    prompt{i+1}=['Ch' int2str(i) ' Ex.'];
end
ChEx = inputdlg(prompt, ...
    'Dataset info (Enter 0 for DIC)', [1 35], DatInf(1:size(prompt,2)));

[UC,~,k] = unique(ChEx);
N = histc(k,1:numel(UC));
if sum(N>1)
ChDup=find(strcmp(UC(N>1),ChEx));
else
    ChDup=[];
end

if ~isempty(ChDup)
    GFPCh=inputdlg({'Which channel is GFP (' int2str(ChDup(1)) ' or ' int2str(ChDup(2))});
    GFPCh=str2double(GFPCh);
end





%read the .lif file and get plane of interest. Indexing with the reader
%object starts with 0, needs to be changed to select a specific series

SeriesI=str2double(ChEx{1});

for j=1:sC
    for i=1:sz
        iPlane=reader.getIndex(i-1,j-1,SeriesI-1)+1; %nSeries is the 0 in getIndex
        image.data(:,:,i,j)=bfGetPlane(reader,iPlane);
    end
end



% Extract the metadata.
hashtable = ser;
keys = arrayfun(@char, hashtable.keySet.toArray, 'UniformOutput', false);
values = cellfun(@(x) hashtable.get(x), keys, 'UniformOutput', false);

% Organize the metadata.
metadata.keys = keys;
metadata.values = values;
metadata.hashtable = hashtable;


zSlices = sz;
%


numChannels=sC;
%
% % Do we have enough color channels?
if numChannels < 3
    f = uifigure;
    str = sprintf('"%s" has less than 3 color channels!', filename);
    uialert(f, str, 'Invalid Image Format');
    return;
end


% Get the image size.
yPixels = sx;
xPixels = sy;

% Get the scale size

xy_s=omeMeta.getPixelsPhysicalSizeX(0).value();
xy_scale=xy_s.doubleValue();
z_s=omeMeta.getPixelsPhysicalSizeZ(0).value();
z_scale=z_s.doubleValue();


% Initialize the image volume information.
image.pixels = [ ...
    xPixels; ...
    yPixels; ...
    zSlices];
image.scale = [ ...
    xy_scale; ...
    xy_scale; ...
    z_scale];

% % Initialize the channel information.

image.channels = cell(numChannels,1);
image.colors = nan(numChannels,3);
image.lasers = nan(numChannels,1);
image.dicChannel=NaN;
%image.emissions = nan(numChannels,1);

for c=1:numChannels
    vEx=str2double(ChEx{c+1});
    if exist('GFPCh','var') && c==GFPCh
        image.channels(c)={'GFP'};
        image.colors(c,:)=[1 1 1];
        image.lasers(c)=vEx;
    elseif vEx==0
        image.channels(c)={'DIC'};
        image.dicChannel=c;
    elseif vEx~=0 && vEx<=410
        image.channels(c)={'mTagBFP'};
        image.colors(c,:)=[0 0 1];
        image.lasers(c)=vEx;
    elseif vEx>=410 && vEx<500
        image.channels(c)={'CyOFP'};
        image.colors(c,:)=[0 1 0];
        image.lasers(c)=vEx;
    elseif vEx>=500 && vEx<590
        image.channels(c)={'TagFRP'};
        image.colors(c,:)=[1 1 1];
        image.lasers(c)=vEx;
    elseif vEx>=590
        image.channels(c)={'mNeptune2.5'};
        image.colors(c,:)=[1 0 0];
        image.lasers(c)=vEx;
    end
end





end
