function o = call_spots_spatial(o)
%% o = o.call_spots_spatial
% 
% This is another probability method for gene calling.
% The difference here, is that an image is built up for each gene and then
% the spots are the local maxima on each gene image. This allows for
% multiple gene matches at each location i.e. overlapping spots. 
%
% o: iss object
% LookupTable: should be returned from call_spots_prob. It
% just gives the the probabilities that each spot score is explained by each
% gene. It saves calculating the probabilities explicitly each time.
%
% produces 
% pxSpotColors(Spot,b,r): intensity of Spot in channel b, round r
% pxSpotGlobalYX(Spot,:): global yx coordinate for each spot
% pxSpotCodeNo(Spot): gene index for each spot
% pxLogProbOverBackground(Spot): log of probability spot can be explained
% by gene relative to probability it can be explained by background.
% pxSpotScore(Spot): pxLogProbOverBackground of best gene match relative to
% second best gene match at that location.
% pxSpotScoreDev(Spot): standard deviation in spot scores across all genes
% at that location.
% pxSpotIntensity(Spot): intensity of the spot. Takes into account
% pxSpotCodeNo. Calculated by get_spot_intensity.
% 
%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end


%% get o.spBledCodes.
if isempty(o.ompBledCodes)
    o = o.get_omp_bled_codes;
end
o.spBledCodes = o.ompBledCodes;

%% Find shape of spots to fit to data
o = o.get_spot_shape;

%% Get genes on each tile
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end
    
AllSpotColors = cell(nTiles,1);
AllLocalYX = cell(nTiles,1);
AllSpotCodeNo = cell(nTiles,1);
AllSpotCoef = cell(nTiles,1);
AllSpotIter = cell(nTiles,1);
AllOriginalTile = cell(nTiles,1);
AllResidueReduction = cell(nTiles,1);

for t=NonemptyTiles
    [LocalYX,SpotColorsInitial] = o.get_spot_colors_all_pixels(t);
    SpotColors = (double(SpotColorsInitial)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    [AllSpotCodeNo{t},AllSpotCoef{t},AllSpotLocalInd,AllSpotIter{t},AllResidueReduction{t}] = ...
        get_spatial_spots(o,SpotColors,LocalYX);
    AllSpotColors{t} = SpotColorsInitial(AllSpotLocalInd,:,:);
    AllLocalYX{t} = LocalYX(AllSpotLocalInd,:);
    AllOriginalTile{t} = ones(size(AllSpotCodeNo{t}))*t;
end
save(fullfile(o.OutputDirectory, 'CallSpotsWorkspace.mat'), 'o','AllSpotColors','AllLocalYX',...
    'AllSpotCodeNo','AllSpotCoef','AllSpotIter','AllOriginalTile','AllResidueReduction');
AllSpotColors = cell2mat(AllSpotColors(NonemptyTiles));
AllLocalYX = cell2mat(AllLocalYX(NonemptyTiles));
AllSpotCodeNo = cell2mat(AllSpotCodeNo(NonemptyTiles));
AllSpotCoef = cell2mat(AllSpotCoef(NonemptyTiles));
AllSpotIter = cell2mat(AllSpotIter(NonemptyTiles));
AllOriginalTile = cell2mat(AllOriginalTile(NonemptyTiles));
AllResidueReduction =  cell2mat(AllResidueReduction(NonemptyTiles));

%% Remove duplicates by keeping only spots detected on their home tile
AllGlobalYX = double(AllLocalYX)+o.TileOrigin(AllOriginalTile,:,o.ReferenceRound);
if o.Graphics
    figure(1350)
    plot(AllGlobalYX(:,2), AllGlobalYX(:,1), '.', 'markersize', 1);
    title('All global coords including duplicates');
    %set(gca, 'YDir', 'reverse');
end
[AllLocalTile, ~] = which_tile(AllGlobalYX, o.TileOrigin(:,:,o.ReferenceRound), o.TileSz);
NotDuplicate = (AllLocalTile==AllOriginalTile);
ndSpotColors = AllSpotColors(NotDuplicate,:,:);
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndSpotCodeNo = AllSpotCodeNo(NotDuplicate,:);
ndSpotCoef= AllSpotCoef(NotDuplicate,:);
ndSpotIter = AllSpotIter(NotDuplicate);
ndOriginalTile = AllOriginalTile(NotDuplicate);
ndResidueReduction = AllResidueReduction(NotDuplicate);

if o.Graphics
    figure(1351)
    plot(ndGlobalYX(:,2), ndGlobalYX(:,1), '.', 'markersize', 1);
    title('Global coords without duplicates');
    drawnow;
    %set(gca, 'YDir', 'reverse');
end
o.spSpotColors = ndSpotColors;
o.spSpotCodeNo = ndSpotCodeNo;
o.spSpotGlobalYX = ndGlobalYX;
o.spLocalTile = ndOriginalTile;
o.spCoefs = ndSpotCoef;
o.spSpotIter = ndSpotIter;
o.spSpotScore = ndResidueReduction;

end
