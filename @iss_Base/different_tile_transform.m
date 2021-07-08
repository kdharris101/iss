function [MyPointCorrectedYX, error, nMatches] = different_tile_transform(o, y0, x0,MyLocalYX, t, t2, r, b)
%% [MyPointCorrectedYX, error, nMatches] = ...
%  o.different_tile_transform(y0, x0,MyLocalYX, t, t2, r, b)
% 
% Using the transformation variables found by the PointCloudRegistration
% algorithm, this tries to match spots which are on tile t on round r but 
% tile t2 on the reference round. Also returns corrected coordinates of
% CenteredMyLocalYX
%
% y0: cell containig the YX location of all spots in all rounds 
% and colour channels for all tiles
% x0: cell containing the YX location of spots in the 
% anchor channel for all tiles
% MyLocalYX: centered local coordinates of spots on tile t on round r
% and tile t2 on the reference (Anchor) round
% MyPointCorrectedYX: YX local coordinates for spots on t2, round r. 
% error: error in the PCR between tile t, reference round and tile t2, round r
% nMatches: matches in the PCR 
%%
y = y0{t,b,r};
x = [vertcat(x0{t2,:})-o.TileCentre,ones(length(vertcat(x0{t2,:})),1)];

if isempty(o.PcDist)
    o.PcDist = inf;
end

% find well isolated points as those whose second neighbor is far
k0 = KDTreeSearcher(y);
[~, d2] = k0.knnsearch(y, 'k', 2);
if isfinite(o.PcDist) && size(y,1) > 1 
    y = y(d2(:,2)>o.PcDist*2,:);
end

%Make kd trees out of these well isolated points
k = KDTreeSearcher(y);
%Transform according to tile in reference round to get transformed local coordinates relative to t2
xM_t2 = x*o.D(:,:,t2,r,b)+o.TileCentre;  
xM_Global = xM_t2 + o.TileOrigin(t2,:,r);     %Add TileOrigin for t2, to get global coordinates
xM_t1 = xM_Global - o.TileOrigin(t,:,r);      %Go from global to local coordinates relative to t1.

[~,Dist] = k.knnsearch(xM_t1);
UseMe = Dist<o.PcDist;               
nMatches = sum(UseMe);
error = sqrt(mean(Dist(UseMe>0).^2));

MyPointCorrectedYX = MyLocalYX*o.D(:,:,t2,r,b)+o.TileOrigin(t2,:,r)-o.TileOrigin(t,:,r);
MyPointCorrectedYX = round(MyPointCorrectedYX+o.TileCentre);

return




