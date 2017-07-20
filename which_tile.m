function [Tile, LocalCoord] = which_tile(Points, TileOrigins, TileSz)
% [Tile, LocalCoord] = which_tile(GlobalCoord, TileOrigins, TileSz)
%
% You have a tiled image that has already been stitched, and the XY coords 
% of the origin of each tile in a global cooridnate system are in
% TileOrigins (size nTiles x 2). All tiles are square of size TileSz
% 
% You have a set of Points in the global coordinate system (nPoints x 2), 
% and you want to know which Tile to find them in (nPoints x 1), and what 
% their LocalCoord in that tile (nPoints x 2)
% 
% sometimes a point is in more than one tile (because overlap) - in this
% case, it will choose the one with the lowest tile index
% 
% NOTE here we assuming coordinates within a tile starting at 0! (matlab usually starts at 1)
% i.e. LocalCoord = GlobalCoord - TileOrigin
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
TileCenters = TileOrigins + TileSz/2 - .5; % .5 because even tile size
nTiles = size(TileCenters,1);
nPoints = size(Points,1);


% l_infinity distance of each point from each tile center (size nPoints x nTiles)
SquareDist = inf(nPoints,nTiles);
for t=1:nTiles
    if all(isfinite(TileCenters(t,:)),2)
        SquareDist(:,t) = max(abs(Points - TileCenters(t,:)),[],2);
    end
end

[MinDist, Tile] = min(SquareDist,[],2);
FitsIn = (MinDist<TileSz/2); % is the point actually in the tile? ASSUME SQUARE TILE!
Tile(~FitsIn) = nan; % if no fit, Tile and LocalCoord is nan
LocalCoord = nan(size(Points));
LocalCoord(FitsIn,:) = Points(FitsIn,:) - TileOrigins(Tile(FitsIn),:);

return
