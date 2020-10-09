function [] = view_filtering(o,r,t)
%%  [] = view_filtering(o,r,t)
%
% Allows you to view filtering on tile t, round r, for all colour channels. 
% You can also see how the radius of filter changes the results using the
% scrollbar.
%
% o: iss object
% r: round
% t: tile
%
% The tiles don't need to have been produced yet but
% o must have the following parameters specified:
% o.InputDirectory, o.FileBase, o.RawFileExtension,
% o.AnchorRound, o.AnchorChannel, o.DapiChannel, o.TileSz  

%% Get paramaters needed for filtering
%Global paramater as plot has multiple uicontrol, makes it easier to share data between them.
global vf_ISSPlotObject         
imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);  %raw data file name for round r 

% construct a Bio-Formats reader with the Memoizer wrapper
bfreader = loci.formats.Memoizer(bfGetReader(), 0);
bfreader.setId(imfile);
% get some basic image metadata
[nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
    get_ome_tilepos(bfreader);
scene = nSeries/nSerieswPos;
bfreader.close();

vf_ISSPlotObject.r = r;
vf_ISSPlotObject.t = t;
vf_ISSPlotObject.nChannels = nChannels;
vf_ISSPlotObject.AnchorRound = o.AnchorRound;
vf_ISSPlotObject.AnchorChannel = o.AnchorChannel;
vf_ISSPlotObject.DapiChannel = o.DapiChannel;
vf_ISSPlotObject.TileSz = o.TileSz;


%% Default Filters parameters - different for Dapi
if strcmpi(o.ExtractR1, 'auto') || strcmpi(o.ExtractR2, 'auto')
    vf_ISSPlotObject.R1 = round(0.5/pixelsize);
    vf_ISSPlotObject.R2 = vf_ISSPlotObject.R1*2;
else
    vf_ISSPlotObject.R1 = o.ExtractR1;
    vf_ISSPlotObject.R2 = o.ExtractR2;
end

if strcmpi(o.DapiR,'auto')
    vf_ISSPlotObject.DapiR = round(8/pixelsize);
else
    vf_ISSPlotObject.DapiR = o.DapiR;
end

%% Default scaling - note 'auto' here is not the same as in extract_and_filter
if r ~=o.AnchorRound
    if strcmpi(o.ExtractScale, 'auto')
        vf_ISSPlotObject.Scale = 2;
    else
        vf_ISSPlotObject.Scale = o.ExtractScale;
    end
elseif r==o.AnchorRound 
    if strcmpi(o.ExtractScaleAnchor, 'auto')
        vf_ISSPlotObject.Scale = 2;
    else
        vf_ISSPlotObject.Scale = o.ExtractScaleAnchor;
    end
end

%% Read in raw tile data
% a new reader per worker
bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
% use the memo file cached before
bfreader.setId(imfile);

bfreader.setSeries(scene*t-1);

vf_ISSPlotObject.RawImages = zeros(o.TileSz,o.TileSz,nChannels);
vf_ISSPlotObject.FilterImages = zeros(o.TileSz,o.TileSz,nChannels);
fprintf('Loading in colour channel   ');
for c=1:nChannels
    fprintf('\b%d',c);
    I = cell(nZstacks,1);
    for z = 1:nZstacks
        iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
        I{z} = bfGetPlane(bfreader, iPlane);
    end
    
    % focus stacking
    I_mod = o.fstack_modified(I(o.FirstZPlane:end));
    vf_ISSPlotObject.RawImages(:,:,c) = I_mod;
end
fprintf('\nFinished loading tile.\n');

%% Plot
viewfilter_plot;
end
