function o = call_spots_pixel(o,LookupTable)
%% o = o.call_spots_pixel(LookupTable)
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
%% Load in images, make images for each gene, find local maxima in images
nCodes = length(o.CharCodes);
rr = o.ReferenceRound;      %Round to base coordinate system on
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end

%Spots on achor round cover whole range of coordinates, same for each tile
AnchorLocalYX = zeros(o.TileSz^2,2);
AnchorLocalYX(:,1) = repelem(1:o.TileSz,1,o.TileSz);
AnchorLocalYX(:,2) = repmat(1:o.TileSz,1,o.TileSz);
nPixels = size(AnchorLocalYX,1);
AnchorLocalYX = [AnchorLocalYX,ones(nPixels,1)];

PeakSpotColors = cell(nCodes,1);
PeakLocalYX = cell(nCodes,1);
PeakLogProbOverBackground = cell(nCodes,1);
Peak2ndBestLogProb = cell(nCodes,1);
PeakScoreDev = cell(nCodes,1);
OriginalTile = cell(nCodes,1);

%Get output file names so don't have data from all tiles in matlab at once
OutputTileNo = cumsum(equal_split(int32(length(NonemptyTiles)),round(length(NonemptyTiles)/o.PixelFileMaxTiles)));
nFiles = length(OutputTileNo);
o.PixelFileNames = cell(nFiles,1);  

for FileIdx=1:nFiles
    if isequal(NonemptyTiles,1:nTiles)
        o.PixelFileNames(FileIdx) =...
            {fullfile(o.OutputDirectory, strcat('ByPixelWorkspace',num2str(OutputTileNo(FileIdx)),'.mat'))};
    else
         o.PixelFileNames(FileIdx) =...
            {fullfile(o.OutputDirectory, strcat('AbridgedByPixelWorkspace',num2str(OutputTileNo(FileIdx)),'.mat'))};
    end
end

FileIdx = 1;
for t=1:length(NonemptyTiles)  
    tile_no = NonemptyTiles(t);
    if exist(o.PixelFileNames{FileIdx}, 'file')
        if ismember(t,OutputTileNo)
            FileIdx=FileIdx+1;
        end
        fprintf('Tile %d already done.\n', tile_no);
        continue;
    end
        
    %Get pixel colors
    [GoodAnchorLocalYX,GoodSpotColors] = o.get_spot_colors(tile_no,AnchorLocalYX,nPixels);
    
    %Get local maxima log probabilities for each gene
    [tPeakLocalYX,tPeakSpotColors,tPeakLogProbOverBackground,...
    tPeak2ndBestLogProb,tPeakScoreDev,tOriginalTile] = ...
    o.detect_peak_genes(LookupTable,GoodSpotColors,GoodAnchorLocalYX,tile_no);
    clearvars GoodSpotColors GoodAnchorLocalYX;
    
    %Keep data for all tiles together
    PeakSpotColors = cellfun( @(x,y) [x;y], PeakSpotColors, tPeakSpotColors, 'UniformOutput', false );
    PeakLocalYX = cellfun( @(x,y) [x;y], PeakLocalYX, tPeakLocalYX, 'UniformOutput', false );
    clearvars tPeakSpotColors tPeakLocalYX;
    PeakLogProbOverBackground = cellfun( @(x,y) [x;y], PeakLogProbOverBackground, tPeakLogProbOverBackground, 'UniformOutput', false );
    Peak2ndBestLogProb = cellfun( @(x,y) [x;y], Peak2ndBestLogProb, tPeak2ndBestLogProb, 'UniformOutput', false );
    PeakScoreDev = cellfun( @(x,y) [x;y], PeakScoreDev, tPeakScoreDev, 'UniformOutput', false );
    OriginalTile = cellfun( @(x,y) [x;y], OriginalTile, tOriginalTile, 'UniformOutput', false );
    clearvars tPeakLogProbOverBackground tPeak2ndBestLogProb tPeakScoreDev tOriginalTile;
    
    if ismember(t,OutputTileNo)
        save(o.PixelFileNames{FileIdx}, 'PeakSpotColors','PeakLocalYX', 'PeakLogProbOverBackground',...
            'Peak2ndBestLogProb','PeakScoreDev','OriginalTile', '-v7.3');
        PeakSpotColors = cell(nCodes,1);
        PeakLocalYX = cell(nCodes,1);
        PeakLogProbOverBackground = cell(nCodes,1);
        Peak2ndBestLogProb = cell(nCodes,1);
        PeakScoreDev = cell(nCodes,1);
        OriginalTile = cell(nCodes,1);
        FileIdx=FileIdx+1;
    end
end

%% Deal with each file one by one
SpotCodeNo = cell(1,1);
SpotColors = cell(1,1);
GlobalYX = cell(1,1);
LogProbOverBackground = cell(1,1);
SecondBestLogProb = cell(1,1);
ScoreDev = cell(1,1);
    
nFiles = length(o.PixelFileNames);
fprintf('\nGetting results from file      ');
for f = 1:nFiles
    if nFiles<10
        fprintf('\b\b\b%d/%d',f,nFiles);
    else
        if f<10; fprintf('\b\b\b\b%d/%d',f,nFiles);
        else; fprintf('\b\b\b\b\b%d/%d',f,nFiles); end
    end

    %Get global coordinates of peaks
    load(cell2mat(o.PixelFileNames(f)));
    PeakGlobalYX = cell(nCodes,1);
    for GeneNo = 1:nCodes
        PeakGlobalYX{GeneNo} = bsxfun(@plus,double(PeakLocalYX{GeneNo}),o.TileOrigin(OriginalTile{GeneNo},:,rr));
    end 
    % Clear memory by removing bad matches, pSpotScore<-5
    clearvars PeakLocalYX
    %Have to filter results initially so don't save too many and have
    %memory issues. Only take gene which are 1st or second best at thier location
    %unless probability relative to background is good.
    QualOK = cellfun(@(x1,x2) x1-x2>=o.pxInitialScoreThresh | x1>o.pxInitialProbThresh,PeakLogProbOverBackground,Peak2ndBestLogProb,'UniformOutput',false);
    PeakSpotColors = cellfun(@(x1,x2) x1(x2,:,:),PeakSpotColors,QualOK,'UniformOutput',false);
    PeakGlobalYX = cellfun(@(x1,x2) x1(x2,:),PeakGlobalYX,QualOK,'UniformOutput',false);
    PeakLogProbOverBackground = cellfun(@(x1,x2) x1(x2),PeakLogProbOverBackground,QualOK,'UniformOutput',false);
    Peak2ndBestLogProb = cellfun(@(x1,x2) x1(x2),Peak2ndBestLogProb,QualOK,'UniformOutput',false);
    PeakScoreDev = cellfun(@(x1,x2) x1(x2),PeakScoreDev,QualOK,'UniformOutput',false);
    OriginalTile = cellfun(@(x1,x2) x1(x2),OriginalTile,QualOK,'UniformOutput',false);
    clearvars QualOK
        
    % Remove duplicates by keeping only spots detected on their home tile
    ndSpotColors = cell(nCodes,1);
    ndGlobalYX = cell(nCodes,1);
    ndLogProbOverBackground = cell(nCodes,1);
    nd2ndBestLogProb = cell(nCodes,1);
    ndScoreDev = cell(nCodes,1);
    
    for GeneNo = 1:nCodes
        if o.Graphics==2
            figure(1001)
            plot(PeakGlobalYX{GeneNo}(:,2), PeakGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
            title('All global coords including duplicates');
            %set(gca, 'YDir', 'reverse');
        end
        
        [AllLocalTile, ~] = which_tile(PeakGlobalYX{GeneNo}, o.TileOrigin(:,:,rr), o.TileSz);
        NotDuplicate = (AllLocalTile==OriginalTile{GeneNo});
        ndSpotColors{GeneNo} = PeakSpotColors{GeneNo}(NotDuplicate,:,:);
        ndGlobalYX{GeneNo} = PeakGlobalYX{GeneNo}(NotDuplicate,:);
        ndLogProbOverBackground{GeneNo} = PeakLogProbOverBackground{GeneNo}(NotDuplicate);
        nd2ndBestLogProb{GeneNo} = Peak2ndBestLogProb{GeneNo}(NotDuplicate);
        ndScoreDev{GeneNo} = PeakScoreDev{GeneNo}(NotDuplicate);
        
        if o.Graphics==2
            figure(1002); clf
            plot(ndGlobalYX{GeneNo}(:,2), ndGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
            title('Global coords without duplicates');
            drawnow;
            %set(gca, 'YDir', 'reverse');
        end
    end
    
    %Free up memory
    clearvars PeakSpotColors PeakGlobalYX PeakLogProbOverBackground Peak2ndBestLogProb PeakScoreDev...
        OriginalTile AllLocalTile NotDuplicate
    
    % Get final results
    nGeneSpots = cell2mat(cellfun(@length,ndScoreDev,'uni',false));    
    for GeneNo = 1:nCodes
        SpotCodeNo{1} = [SpotCodeNo{1};ones(nGeneSpots(GeneNo),1)*GeneNo];
        SpotColors{1} = [SpotColors{1};ndSpotColors{GeneNo}];
        GlobalYX{1} = [GlobalYX{1};ndGlobalYX{GeneNo}];
        LogProbOverBackground{1} = [LogProbOverBackground{1};ndLogProbOverBackground{GeneNo}];
        SecondBestLogProb{1} = [SecondBestLogProb{1};nd2ndBestLogProb{GeneNo}];
        ScoreDev{1} = [ScoreDev{1};ndScoreDev{GeneNo}];
    end
    clearvars ndSpotColors ndGlobalYX ndLogProbOverBackground nd2ndBestLogProb ndScoreDev
end
fprintf('\n');

%% Add results to iss object

o.pxSpotColors = cell2mat(SpotColors);
o.pxSpotCodeNo = cell2mat(SpotCodeNo);
o.pxSpotGlobalYX = cell2mat(GlobalYX);
o.pxLogProbOverBackground = cell2mat(LogProbOverBackground);
o.pxSpotScore = o.pxLogProbOverBackground-cell2mat(SecondBestLogProb);
o.pxSpotScoreDev = cell2mat(ScoreDev);        
o.pxSpotIntensity = o.get_spot_intensity(o.pxSpotCodeNo,o.pxSpotColors);
end


