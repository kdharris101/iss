function [RoundTile,RoundYX] = get_SpotTileEachRound(o,SpotGlobalYX,LocalTile)
%% decide which tile to read each spot off in each round. 
% They are read of home tile if possible (always possible in ref round)
% in other rounds might have to be a NWSE neighbor - but never a diagonal
% neighbor
% SpotGlobalYX: anchor round coordinates of spots
% LocalTile: tile each spot was found on
% ndRoundTile(s,r) stores appropriate tile for spot s on round r
% ndRoundYX(s,:,r) stores YX coord on this tile
nSpots = size(LocalTile,1);
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
RoundTile = nan(nSpots,max(o.UseRounds));
RoundYX = nan(nSpots,2,max(o.UseRounds));
PossNeighbs = [-1 -nY 1 nY 0];   %NWSE then same tile - same will have priority by being last

for r=o.UseRounds
    if length(LocalTile)>10 && size(SpotGlobalYX,1)>1000
        fprintf('Finding appropriate tiles for round %d\n', r);
    end
    
    for n = PossNeighbs
        % find origins of each tile's neighbor, NaN if not there
        NeighbTile = (1:nTiles)+n;
        NeighbOK = (NeighbTile>=1 & NeighbTile<=nTiles);
        NeighbOrigins = nan(nTiles,2);
        NeighbOrigins(NeighbOK,:) = round(o.TileOrigin(NeighbTile(NeighbOK),:,r));
        
        % now for each spot see if it is inside neighbor's tile area
        SpotsNeighbOrigin = NeighbOrigins(LocalTile,:);
        SpotsInNeighbTile = all(SpotGlobalYX>=SpotsNeighbOrigin+1+o.ExpectedAberration...
            & SpotGlobalYX<=SpotsNeighbOrigin+o.TileSz-o.ExpectedAberration, 2);
        
        % for those that were in set this to be its neighbor
        RoundTile(SpotsInNeighbTile,r) = NeighbTile(LocalTile(SpotsInNeighbTile));    
    end
    
    % compute YX coord
    HasTile = isfinite(RoundTile(:,r));
    RoundYX(HasTile,:,r) = SpotGlobalYX(HasTile,:) - ...
        round(o.TileOrigin(RoundTile(HasTile,r),:,r));
    
end
end

