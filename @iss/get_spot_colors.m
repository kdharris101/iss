function [GoodAnchorLocalYX,GoodSpotColors] = get_spot_colors(o,t,AnchorLocalYX,nPixels)
%%
% This  returns the spot colors across all pixels in each
% round/channel
% Input
% o - iss object
% t - tile of interest
% AnchorLocalYX(:,1:2) - YX coordinates covering whole of tile.
% AnchorLocalYX(:,3) = 1
% nPixels - just length of AnchorLocalYXl
% Output
% GoodAnchorLocalYX - YX coordinates in anchor image that when shifted stay
% within range for all rounds and channels
% GoodSpotColors - GoodSpotColors(:,b,r) is the spot intensity in channel
% b, round r at anchor locations given by GoodAnchorLocalYX

%%
fprintf('Tile %d: Getting spot colours\n', t);

TileSpotColors = nan(nPixels, o.nBP, o.nRounds);

for r = o.UseRounds
    FileName = o.TileFiles{r,t};
    TifObj = Tiff(FileName);
    for b = o.UseChannels
        TifObj.setDirectory(o.FirstBaseChannel + b - 1);
        BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
        %Find shifted shifted coordinates on new round/channel and
        %corresponding intensities
        MyPointCorrectedYX = o.A(b)*(AnchorLocalYX*o.D(:,:,t,r));
        MyPointCorrectedYX = round(MyPointCorrectedYX);
        TileSpotColors(:,b,r) = IndexArrayNan(BaseIm, MyPointCorrectedYX');
    end
end

%Only keep spots that are in all rounds/channels
SpotColorsToUse = TileSpotColors(:,o.UseChannels,o.UseRounds);
Good = all(isfinite(SpotColorsToUse(:,:)),2);
GoodAnchorLocalYX = AnchorLocalYX(Good,1:2);
GoodSpotColors = TileSpotColors(Good,:,:);