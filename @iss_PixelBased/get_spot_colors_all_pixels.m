function [GoodAnchorLocalYX,GoodSpotColors] = get_spot_colors_all_pixels(o,t,AnchorLocalYX)
%%
% This  returns the spot colors across all pixels in each
% round/channel
% Input
% o - iss object
% t - tile of interest
% AnchorLocalYX - yx position of coordinates of interest on tile t. 
% Output
% GoodAnchorLocalYX - YX coordinates in anchor image that when shifted stay
% within range for all rounds and channels
% GoodSpotColors - GoodSpotColors(:,b,r) is the spot intensity in channel
% b, round r at anchor locations given by GoodAnchorLocalYX

%%
if nargin<3 || isempty(AnchorLocalYX)
    fprintf('Tile %d: Getting spot colours\n', t);
    %Spots on achor round cover whole range of coordinates, same for each tile
    AnchorLocalYX = zeros(o.TileSz^2,2);
    AnchorLocalYX(:,1) = repelem(1:o.TileSz,1,o.TileSz);
    AnchorLocalYX(:,2) = repmat(1:o.TileSz,1,o.TileSz);
end
nPixels = size(AnchorLocalYX,1);
AnchorLocalYX = [AnchorLocalYX,ones(nPixels,1)];
    
CenteredAnchorLocalYX = AnchorLocalYX - [o.TileCentre,0];
TileSpotColors = nan(nPixels, o.nBP, o.nRounds);

%Maybe need to include case where tile changes with rounds!!
for r = o.UseRounds
    FileName = o.TileFiles{r,t};
    TifObj = Tiff(FileName);
    for b = o.UseChannels
        TifObj.setDirectory(o.FirstBaseChannel + b - 1);
        BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
        %Find shifted shifted coordinates on new round/channel and
        %corresponding intensities
        MyPointCorrectedYX = CenteredAnchorLocalYX*o.D(:,:,t,r,b)+o.TileCentre;
        MyPointCorrectedYX = round(MyPointCorrectedYX);
        TileSpotColors(:,b,r) = IndexArrayNan(BaseIm, MyPointCorrectedYX');
    end
end

%Only keep spots that are in all rounds/channels
SpotColorsToUse = TileSpotColors(:,o.UseChannels,o.UseRounds);
Good = all(isfinite(SpotColorsToUse(:,:)),2);
GoodAnchorLocalYX = int32(AnchorLocalYX(Good,1:2));
GoodSpotColors = int32(TileSpotColors(Good,:,:));