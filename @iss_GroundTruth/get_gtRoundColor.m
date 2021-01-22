function o = get_gtRoundColor(o,Method)
%% o = o.get_gtRoundColor(Method);
% Find intensity of all spots found by 'Method' in the gt round
if ~ismember({Method},o.CallMethods)
    error('Method invalid, must be member of o.CallMethods');
end
pf = o.CallMethodPrefix(Method);

%Basic Info
UseRoundsOrig = o.UseRounds;
o.UseRounds = o.gtRounds;
SpotGlobalYX = o.([pf,'SpotGlobalYX']);
LocalTile = o.([pf,'LocalTile']);
SpotLocalYX = SpotGlobalYX-o.TileOrigin(LocalTile,:,o.ReferenceRound);

[RoundTile,~] = o.get_SpotTileEachRound(SpotGlobalYX,LocalTile);
% loop through all tiles, finding spot colors
[SpotColors,~] = get_spot_colors(o,SpotLocalYX,LocalTile,...
    RoundTile,o.gtRawLocalYX);
o.([pf,'_gtColor']) = SpotColors;
o.UseRounds = UseRoundsOrig;
end

