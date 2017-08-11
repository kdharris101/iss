function [o, CellMap, DapiBoundaries] = segment_dapi(o, Dapi0)
% [0 CellMap BoundaryImage] = o.segment_dapi(DapiIm)
%
% segments a DAPI image and assigns each pixel to a cell. Input is
% optional, otherwise it will load from o.BigDapiFile
%
% only works within region outlined by o.CellCallRegionYX
%
% Output CellMap is same size as input DapiIm, with integer entries for
% each pixel, saying which cell it belongs to. (Zero if unassigned)
%
% also saved to o.CellMapFile


%% find Cell Calling Region
y0 = min(o.CellCallRegionYX(:,1));
x0 = min(o.CellCallRegionYX(:,2));
y1 = max(o.CellCallRegionYX(:,1));
x1 = max(o.CellCallRegionYX(:,2));

Mask = poly2mask(o.CellCallRegionYX(:,2)-x0+1, o.CellCallRegionYX(:,1)-y0+1, y1-y0+1, x1-x0+1);
if nargin<2
    Dapi = imread(o.BigDapiFile, 'PixelRegion', {[y0 y1], [x0 x1]}).*uint16(Mask);
else
    Dapi = Dapi0(y0:y1, x0:x1) .*uint16(Mask);
end

%%
Dapi = imadjust(Dapi); % contrast enhancement
ImSz = size(Dapi);
Debug = 1;
%% threshold the map
ThreshVal = prctile(Dapi(Mask), o.DapiThresh);

bwDapi = imerode(Dapi>ThreshVal, strel('disk', 2));

if Debug
    figure(300)
    subplot(2,1,1);
    imagesc(Dapi); 
    subplot(2,1,2);
    imagesc(bwDapi);
    colormap bone
    fprintf('Threshold = %f\n', ThreshVal);
    
end
%% find local maxima 
dist = bwdist(~bwDapi);
dist0 = dist;
dist0(dist<o.DapiMinSize)=0;
ddist = imdilate(dist0, strel('disk', o.DapiMinSep));
%clear dist 
impim = imimposemin(-dist0, imregionalmax(ddist));
clear dist0
if Debug
    figure(301);
    subplot(2,1,1)
    imagesc(dist);
    subplot(2,1,2)
    imagesc(impim);
end
%% segment
% remove pixels at watershed boundaries
bwDapi0 = bwDapi;
bwDapi0(watershed(impim)==0)=0;

% assign all pixels a label
labels = uint32(bwlabel(bwDapi0));
[d, idx] = bwdist(bwDapi0);

% now expand the regions by a margin
CellMap0 = zeros(ImSz, 'uint32');
Expansions = (d<o.DapiMargin);
CellMap0(Expansions) = labels(idx(Expansions));

% get rid of cells that are too small
rProps0 = regionprops(CellMap0); % returns XY coordinate and area
BigEnough = [rProps0.Area]>o.MinCellArea;
NewNumber = zeros(length(rProps0),1);
NewNumber(~BigEnough) = 0;
NewNumber(BigEnough) = 1:sum(BigEnough);
CellMap = CellMap0;
CellMap(CellMap0>0) = NewNumber(CellMap0(CellMap0>0));

if Debug
    figure(302)
    image(label2rgb(CellMap, 'jet', 'w', 'shuffle'));

end

%%
% CellYX = fliplr(vertcat(rProps(BigEnough).Centroid)); % because XY

%% make image with boundaries
Boundaries = (CellMap ~= imdilate(CellMap,strel('disk', 1)));
DapiBoundaries = Dapi;

OddPix = mod((1:size(Dapi,2)) + (1:size(Dapi,1))', 2);
DapiBoundaries(Boundaries & OddPix) = .3 * max(Dapi(:));
DapiBoundaries(Boundaries & ~OddPix) = 0;

imwrite(DapiBoundaries, [o.OutputDirectory '\background_boundaries.tif']);

o.CellMapFile = fullfile(o.OutputDirectory, 'CellMap.mat');
save(o.CellMapFile, 'CellMap', 'DapiBoundaries', 'y0', 'y1', 'x0', 'x1');


end

