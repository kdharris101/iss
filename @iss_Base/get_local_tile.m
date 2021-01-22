function t = get_local_tile(o,yx)
%% t = o.get_local_tile(yx)
%returns tile that each in set of global yx coordinates was found on.
%o: iss object
%yx(s,:): global yx coordinates of s locations of interests.
%
nSpots = size(yx,1);
nTiles = size(o.TileOrigin,1);
YXDist2Tiles = repelem(yx,nTiles,1)-...
    repmat(o.TileOrigin(:,:,o.ReferenceRound),nSpots,1);
%If tile in wrong direction then set distance to infinity so won't be min
%distance.
YXDist2Tiles(YXDist2Tiles(:,1)<0 | YXDist2Tiles(:,2)<0,1) = inf;
Dist2Tiles = reshape(vecnorm(YXDist2Tiles,2,2),[nTiles,nSpots])';
[~,t] = min(Dist2Tiles,[],2);
end

