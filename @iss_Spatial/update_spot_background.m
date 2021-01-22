function SpotColors = ...
    update_spot_background(o,SpotColors,BackgroundCodeNo)
%% SpotColors = ...
%   update_spot_background(o,SpotColors,BackgroundCodeNo)
% Update Spot Colors by finding how background can best represent each
% pixel and then subtracting this.
% o: iss object
% SpotColors: SpotColors(s,b,r) is the z-scored intensity in channel b,
% round r at LocalYX(s,:).
% BackgroundCodeNo: Index of background vectors in o.spBledCodes. Must
% be orthogonal.

pf = o.CallMethodPrefix('Spatial');
nSpots = size(SpotColors,1);
DotProduct = SpotColors(:,:)*o.([pf,'BledCodes'])(BackgroundCodeNo,:)';
ToRemove = reshape(DotProduct*o.([pf,'BledCodes'])(BackgroundCodeNo,:),nSpots,o.nBP,o.nRounds);
SpotColors = SpotColors - ToRemove;

end
