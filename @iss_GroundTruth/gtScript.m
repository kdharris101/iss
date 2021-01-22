%% Stuff to check
%o.nExtraRounds: 1 if no anchor, 2 if using anchor. (I.e. number of
%non-imaging rounds)
%o.FileBase(o.gtRound) needs to contain Ground truth raw data.
%o.TileFiles needs to include o.gtRound
%o.AutoThresh needs to include o.gtRound

%Set o.gtGeneNo(b) to be the gene number of the gene present in channel b of
%the ground truth round. 0 means no gene is present. 
o.gtGeneNo = zeros(o.nRounds+o.nExtraRounds,o.nBP);
o.gtGeneNo(8,4) = 22;
o.gtGeneNo(9,2) = 50;
o.gtGeneNo(9,4) = 31;
o.gtGeneNo(9,6) = 64;
o.gtGeneNo(9,7) = 19;

o.gtReferenceChannel = zeros(o.nRounds+o.nExtraRounds,1);
o.gtReferenceChannel(8) = o.AnchorChannel;
o.gtReferenceChannel(9) = 1;
o.gtRounds = find(o.gtReferenceChannel>0)';
o.gtAnchorChannel = zeros(o.nRounds+o.nExtraRounds,1);
o.gtAnchorChannel(8) = o.gtReferenceChannel(8);
o.gtAnchorChannel(9) = o.gtReferenceChannel(9);

UseRoundsOrig = o.UseRounds;
if max(UseRoundsOrig)>o.nRounds
    UseRoundsOrig = 1:o.nRounds;
    warning('o.UseRounds contains round outside 1:o.nRounds so setting o.UseRounds=1:o.nRounds');
end
o.UseRounds = o.gtRounds;

%% Get ground truth Round spots
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
o.TileCentre = 0.5*[o.TileSz+1,o.TileSz+1];
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end

o.gtRawLocalYX = cell(nTiles,o.nBP, o.nRounds+o.nExtraRounds);
o.gtRawIsolated = cell(nTiles,o.nBP, o.nRounds+o.nExtraRounds);

for t=NonemptyTiles
    fprintf('Getting tile %d ground truth spots\n', t);
    for r=o.UseRounds
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        for b=o.UseChannels
            TifObj.setDirectory(o.FirstBaseChannel + b - 1);
            ReferenceIm = int32(TifObj.read())-o.TilePixelValueShift;
            if o.AutoThresh(t,b,r)==0
                o.AutoThresh(t,b,r) = median(abs(ReferenceIm(:)))*o.AutoThreshMultiplier;
            end
            if o.gtGeneNo(r,b)~=0
                %If have ground truth data, must collect all spots with
                %intensity above o.gtColorFalsePositiveThresh.
                o.AutoThresh(t,b,r) = min(o.AutoThresh(t,b,r),o.gtColorFalsePositiveThresh);
            end
            if o.SmoothSize
                SE = fspecial('disk', o.SmoothSize);
                ReferenceImSm = imfilter(ReferenceIm ,SE);
            else
                ReferenceImSm = ReferenceIm;
            end
            [o.gtRawLocalYX{t,b,r}, o.gtRawIsolated{t,b,r}] = ...
                o.detect_spots(ReferenceImSm,t,b,r);
        end
    end
end


%% Align to anchor round
o.gtRegInfo.Scores = zeros(nTiles,max(o.UseRounds));
o.gtRegInfo.ChangedSearch = 0;

for t=NonemptyTiles
    for r = o.UseRounds
        if r == o.AnchorRound
            o.D0(t,:,r) = [0,0];
        else
            [o.D0(t,:,r), o.gtRegInfo.Scores(t,r),tChangedSearch] = ...
                o.get_initial_shift2(o.gtRawLocalYX{t,o.gtReferenceChannel(r),r},...
                vertcat(o.RawLocalYX{t,:}), o.FindSpotsSearch{3},'FindSpots');
            o.gtRegInfo.ChangedSearch = o.gtRegInfo.ChangedSearch+tChangedSearch;
            
            fprintf('Tile %d, shift from reference round (%d) to round %d: [%d %d], score %d \n',...
                t, o.ReferenceRound,r, o.D0(t,:,r),o.gtRegInfo.Scores(t));
            
            %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
            if t == 3 || (mod(o.gtRegInfo.ChangedSearch,2) == 0) && (o.gtRegInfo.ChangedSearch>0)
                o = o.GetNewSearchRange_FindSpots(t,3);
            end
        end
    end
end

if o.PointCloudMethod==1
    o = o.PointCloudRegister6_GT(o.gtRawLocalYX, o.RawLocalYX, nTiles);
elseif o.PointCloudMethod==2
    o = o.PointCloudRegister2_GT(o.gtRawLocalYX, o.RawLocalYX, nTiles);
else
    error('Require o.PointCloudMethod to be 1 or 2');
end


%% Save background images
o.gtBigImFiles = cell(o.nRounds+o.nExtraRounds,o.nBP);
for r=o.UseRounds
    %Compute approx new shifts from D matrices
    YXgtRoundTileShifts = permute(squeeze(o.D(3,:,:,r,o.gtReferenceChannel(r))),[2,1,3]);
    o.TileOrigin(:,:,r) =  o.TileOrigin(:,:,o.ReferenceRound) - YXgtRoundTileShifts;
    
    AnchorOrigin = round(o.TileOrigin(:,:,r));
    MaxTileLoc = max(AnchorOrigin);
    NegOriginShift = abs(min(min(AnchorOrigin(:)),1)-1);        %To correct for negative origin - remove later
    %To correct for spots being outside image bounds i.e. zero pad:
    MissedSpotsShift = abs(min(min(ceil((MaxTileLoc + o.TileSz) - max(o.pxSpotGlobalYX))),0));
    gtBigIm = zeros(ceil([(MaxTileLoc + o.TileSz + NegOriginShift+MissedSpotsShift),o.nBP]), 'uint16');
    
    for t=NonemptyTiles
        MyOrigin = AnchorOrigin(t,:)+NegOriginShift;
        if mod(t,10)==0; fprintf('Loading tile %d ground truth image\n', t); end
        if ~isfinite(MyOrigin(1)); continue; end
        for b=1:o.nBP
            if o.gtGeneNo(r,b)==0; continue; end
            LocalgtIm = imread(o.TileFiles{r,t}, b);
            gtBigIm(floor(MyOrigin(1))+(1:o.TileSz), ...
                floor(MyOrigin(2))+(1:o.TileSz), b) ...
                = imresize(LocalgtIm, 1);
        end
    end
    
    for b=1:o.nBP
        if o.gtGeneNo(r,b)==0; continue; end
        o.gtBigImFiles{r,b} = fullfile(o.OutputDirectory, [o.GeneNames{o.gtGeneNo(r,b)}, '_image.tif']);
        imwrite(gtBigIm(1+NegOriginShift:end,1+NegOriginShift:end,b), o.gtBigImFiles{r,b});
    end
end

%% Find intensity of all spots found by different methods in the ground truth round
o = o.get_gtRoundColor('Prob');
o = o.get_gtRoundColor('Pixel');
o = o.get_gtRoundColor('OMP');
o = o.get_gtRoundColor('Spatial');
%% Take ground truth local maxima as the spots - find SpotColors for these
%For each ground truth channel containing gene, get spot colors.

%Need AllBaseLocalYX in workspace to do different_tile_transform.
load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
    
o.gtSpotGlobalYX = cell(o.nRounds+o.nExtraRounds,o.nBP);
o.gtSpotColors = cell(o.nRounds+o.nExtraRounds,o.nBP);
o.gtLocalTile = cell(o.nRounds+o.nExtraRounds,o.nBP);
o.gt_gtColor = cell(o.nRounds+o.nExtraRounds,o.nBP);
for r=o.UseRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        %Get Anchor round coordinates of ground truth Peaks
        nAll = sum(sum(cellfun(@numel, o.gtRawIsolated(:,b,r))));
        AllGlobalYX = zeros(nAll,2);
        AllLocalYX = zeros(nAll,2);
        AllIsolated = zeros(nAll,1);
        OriginalTile = zeros(nAll,1);
        
        ind = 1;
        for t=NonemptyTiles
            nMySpots = sum(sum(cellfun(@numel, o.gtRawIsolated(t,b,r))));
            AllLocalYX(ind:ind+nMySpots-1,:) = (o.gtRawLocalYX{t,b,r}-...
                o.TileCentre-o.D(3,:,t,r,b))/...
                o.D(1:2,:,t,r,b)+o.TileCentre;
            AllGlobalYX(ind:ind+nMySpots-1,:) = AllLocalYX(ind:ind+nMySpots-1,:)+...
                o.TileOrigin(t,:,o.ReferenceRound);
            OriginalTile(ind:ind+nMySpots-1) = t;
            ind = ind+nMySpots;
        end
        
        %get intensity in all channels of ground truth round.
        o.UseRounds = o.gtRounds;
        [RoundTile,~] = o.get_SpotTileEachRound(AllGlobalYX,OriginalTile);
        RoundTile(:,r) = OriginalTile;      %As we know the tile in round r for certain.
        [All_gtColor,~] = get_spot_colors(o,AllLocalYX,OriginalTile,RoundTile,...
            o.gtRawLocalYX);
        All_gtColor = All_gtColor(:,:,:);
        
        %Remove duplciates
        [AllLocalTile, ~] = which_tile(AllGlobalYX, o.TileOrigin(:,:,o.ReferenceRound), o.TileSz);
        NotDuplicate = (AllLocalTile==OriginalTile);
        ndGlobalYX = AllGlobalYX(NotDuplicate,:);
        ndLocalYX = AllLocalYX(NotDuplicate,:);
        ndLocalTile = AllLocalTile(NotDuplicate,:);
        nd_gtColor = All_gtColor(NotDuplicate,:,:);
        
        %get intensity in all channels of imaging rounds.
        o.UseRounds = UseRoundsOrig;
        [ndRoundTile,~] = o.get_SpotTileEachRound(ndGlobalYX,ndLocalTile);
        [ndSpotColors,~] = get_spot_colors(o,ndLocalYX,ndLocalTile,...
            ndRoundTile,AllBaseLocalYX,o.nMatches,o.Error);
        
        ndSpotColorsToUse = ndSpotColors(:,o.UseChannels,o.UseRounds);
        Good = all(isfinite(ndSpotColorsToUse(:,:)),2);
        o.gtSpotGlobalYX{r,b} = ndGlobalYX(Good,:);
        o.gtSpotColors{r,b} = ndSpotColors(Good,:,:);
        o.gtLocalTile{r,b} = ndLocalTile(Good);
        o.gt_gtColor{r,b} = nd_gtColor(Good,:,:);
    end
end

%% Get PixelBased probabilities for ground truth spots
%Need LookupTable to get PixelBased probabilities.
load(fullfile(o.OutputDirectory, ['LookupTable',num2str(o.ProbMethod),'.mat']));

o.gt_pxSpotGlobalYX = cell(sum(o.gtGeneNo(:)>0),1);
o.gt_pxSpotCodeNo = cell(sum(o.gtGeneNo(:)>0),1);
o.gt_pxLogProbOverBackground = cell(sum(o.gtGeneNo(:)>0),1);
o.gt_pxSpotScoreDev = cell(sum(o.gtGeneNo(:)>0),1);
o.gt_pxSpotScore = cell(sum(o.gtGeneNo(:)>0),1);
o.gt_pxSpotIntensity = cell(sum(o.gtGeneNo(:)>0),1);

i=1;
o.UseRounds = UseRoundsOrig;
for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        LogProbOverBackground = o.get_LogProbOverBackground(o.gtSpotColors{r,b},LookupTable);
        o.gt_pxLogProbOverBackground{i} = LogProbOverBackground(:,o.gtGeneNo(r,b));
        o.gt_pxSpotScoreDev{i} = std(LogProbOverBackground,[],2);
        %Find 2nd best gene so can give score relative to it
        nPeaks = size(o.gt_pxSpotScoreDev{i},1);
        [~,gInd] = max(LogProbOverBackground,[],2);
        LogProbOverBackground(sub2ind(size(LogProbOverBackground),(1:nPeaks)',gInd)) = -inf;
        SecondBestLogProb = max(LogProbOverBackground,[],2);
        o.gt_pxSpotScore{i} =  o.gt_pxLogProbOverBackground{i}-SecondBestLogProb;
        o.gt_pxSpotCodeNo{i} = ones(nPeaks,1)*o.gtGeneNo(r,b);
        o.gt_pxSpotIntensity{i} = ...
            o.get_spot_intensity(o.gt_pxSpotCodeNo{i},o.gtSpotColors{r,b});
        o.gt_pxSpotGlobalYX{i} = o.gtSpotGlobalYX{r,b};
        i=i+1;
    end
end
o.gt_pxSpotGlobalYX = cell2mat(o.gt_pxSpotGlobalYX);
o.gt_pxSpotCodeNo = cell2mat(o.gt_pxSpotCodeNo);
o.gt_pxLogProbOverBackground = cell2mat(o.gt_pxLogProbOverBackground);
o.gt_pxSpotScoreDev = cell2mat(o.gt_pxSpotScoreDev);
o.gt_pxSpotScore = cell2mat(o.gt_pxSpotScore);
o.gt_pxSpotIntensity = cell2mat(o.gt_pxSpotIntensity);

%% Get OMP results for ground truth spots
o.gt_ompCoefs = cell(sum(o.gtGeneNo(:)>0),1);
o.gt_ompSpotScore = cell(sum(o.gtGeneNo(:)>0),1);

i=1;
for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        nPeaks = size(o.gtSpotColors{r,b},1);
        SpotColors = (o.gtSpotColors{r,b}-o.z_scoreSHIFT)./o.z_scoreSCALE;
        o.gt_ompCoefs{i} = o.get_omp_coefs(SpotColors);
        o.gt_ompSpotScore{i} = o.get_omp_score(SpotColors,o.gt_ompCoefs{i},...
            ones(nPeaks,1)*o.gtGeneNo(r,b));
        i=i+1;
    end
end
o.gt_ompCoefs = cell2mat(o.gt_ompCoefs);
o.gt_ompSpotScore = cell2mat(o.gt_ompSpotScore);

%% Get true positive and false positive stats for prob and pixel based methods
o = get_pf_gtTruePositiveSets(o,'Prob');
o = get_pf_gtTruePositiveSets(o,'Pixel');
o = get_pf_gtTruePositiveSets(o,'OMP');
o = get_pf_gtTruePositiveSets(o,'Spatial');