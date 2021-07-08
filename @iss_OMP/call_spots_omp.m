function o = call_spots_omp(o)
%% o = o.call_spots_omp
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
%% 
if isempty(o.ompBledCodes)
    o = o.get_omp_bled_codes;
end
%% Load in images, run OMP on each pixel, find local maxima in omp coefficient 
% images for each gene.
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end

nCodes = length(o.CharCodes);
PeakSpotColors = cell(nCodes,1);
PeakLocalYX = cell(nCodes,1);
PeakCoefs = cell(nCodes,1);
PeakNeighbNonZeros = cell(nCodes,1);
OriginalTile = cell(nCodes,1);
fprintf('Doing OMP\n');
for t=1:length(NonemptyTiles)  
    tile_no = NonemptyTiles(t);
    %Get pixel colors
    [AllAnchorLocalYX,AllSpotColors] = o.get_spot_colors_all_pixels(tile_no);
    AllSpotColors = (double(AllSpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    coefs = o.get_omp_coefs(AllSpotColors);
    
    %Get local maxima in coefs for each gene
    [tPeakLocalYX,tPeakSpotColors,tPeakCoefs,tPeakNeighbNonZeros,tOriginalTile] = ...
    o.detect_peak_genes_omp(coefs,AllSpotColors,AllAnchorLocalYX,tile_no);
    clearvars GoodSpotColors GoodAnchorLocalYX;
    
    %Keep data for all tiles together
    PeakSpotColors = cellfun( @(x,y) [x;y], PeakSpotColors, tPeakSpotColors, 'UniformOutput', false );
    PeakLocalYX = cellfun( @(x,y) [x;y], PeakLocalYX, tPeakLocalYX, 'UniformOutput', false );
    clearvars tPeakSpotColors tPeakLocalYX;
    PeakCoefs = cellfun( @(x,y) [x;y], PeakCoefs, tPeakCoefs, 'UniformOutput', false );
    PeakNeighbNonZeros = cellfun( @(x,y) [x;y], PeakNeighbNonZeros, tPeakNeighbNonZeros, 'UniformOutput', false );
    OriginalTile = cellfun( @(x,y) [x;y], OriginalTile, tOriginalTile, 'UniformOutput', false );
    clearvars tPeakCoefs tPeakNeighbourhoodNonZeros tOriginalTile;    
end


%% Remove duplicates by keeping only spots detected on their home tile
PeakGlobalYX = cell(nCodes,1);
for GeneNo = 1:nCodes
    PeakGlobalYX{GeneNo} = bsxfun(@plus,double(PeakLocalYX{GeneNo}),o.TileOrigin(OriginalTile{GeneNo},:,o.ReferenceRound));
end 

ndSpotColors = cell(nCodes,1);
ndGlobalYX = cell(nCodes,1);
ndCoefs = cell(nCodes,1);
ndNeighbNonZeros = cell(nCodes,1);
ndOriginalTile = cell(nCodes,1);

for GeneNo = 1:nCodes
    if o.Graphics==2
        figure(1001)
        plot(PeakGlobalYX{GeneNo}(:,2), PeakGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
        title('All global coords including duplicates');
        %set(gca, 'YDir', 'reverse');
    end

    [PeakLocalTile, ~] = which_tile(PeakGlobalYX{GeneNo}, o.TileOrigin(:,:,o.ReferenceRound), o.TileSz);
    NotDuplicate = (PeakLocalTile==OriginalTile{GeneNo});
    ndSpotColors{GeneNo} = PeakSpotColors{GeneNo}(NotDuplicate,:,:);
    ndGlobalYX{GeneNo} = PeakGlobalYX{GeneNo}(NotDuplicate,:);
    ndCoefs{GeneNo} = PeakCoefs{GeneNo}(NotDuplicate,:);
    ndNeighbNonZeros{GeneNo} = PeakNeighbNonZeros{GeneNo}(NotDuplicate);
    ndOriginalTile{GeneNo} = OriginalTile{GeneNo}(NotDuplicate);

    if o.Graphics==2
        figure(1002); clf
        plot(ndGlobalYX{GeneNo}(:,2), ndGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
        title('Global coords without duplicates');
        drawnow;
        %set(gca, 'YDir', 'reverse');
    end
end

%Free up memory
clearvars PeakSpotColors PeakGlobalYX PeakCoefs PeakNeighbNonZeros ...
    OriginalTile AllLocalTile NotDuplicate

%% Get final results
SpotCodeNo = cell(1,1);
SpotColors = cell(1,1);
GlobalYX = cell(1,1);
Coefs = cell(1,1);
NeighbNonZeros = cell(1,1);
GoodOriginalTile = cell(1,1);
nGeneSpots = cell2mat(cellfun(@length,ndNeighbNonZeros,'uni',false));    
for GeneNo = 1:nCodes
    SpotCodeNo{1} = [SpotCodeNo{1};ones(nGeneSpots(GeneNo),1)*GeneNo];
    SpotColors{1} = [SpotColors{1};ndSpotColors{GeneNo}];
    GlobalYX{1} = [GlobalYX{1};ndGlobalYX{GeneNo}];
    Coefs{1} = [Coefs{1};ndCoefs{GeneNo}];
    NeighbNonZeros{1} = [NeighbNonZeros{1};ndNeighbNonZeros{GeneNo}];
    GoodOriginalTile{1} = [GoodOriginalTile{1};ndOriginalTile{GeneNo}];
end
clearvars ndSpotColors ndGlobalYX ndCoefs ndNeighbNonZeros ndOriginalTile


%% Add results to iss object

o.ompSpotColors = cell2mat(SpotColors);
o.ompSpotColors = o.ompSpotColors.*o.z_scoreSCALE + o.z_scoreSHIFT;
o.ompSpotCodeNo = cell2mat(SpotCodeNo);
o.ompSpotGlobalYX = cell2mat(GlobalYX);
o.ompLocalTile = cell2mat(GoodOriginalTile);  
o.ompCoefs = cell2mat(Coefs);
o.ompNeighbNonZeros = cell2mat(NeighbNonZeros);
o.ompSpotIntensity = o.get_spot_intensity(o.ompSpotCodeNo,o.ompSpotColors);
o.ompSpotScore = o.get_omp_score;
end

