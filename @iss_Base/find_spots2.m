function o = find_spots2(o, AllBaseLocalYX, SkipRegistration)        
%% o = o.find_spots2(AllBaseLocalYX, SkipRegistration);
%
% finds spots in all tiles using the reference channel, removes
% duplicates in overlap regions and returns nSpots x 2 array o.SpotGlobalYX of
% coordinates in global frame
% 
% Looks up colors from apporpriate files and makes nSpots x nBP x nRounds
% array o.cSpotColors
%
% o.Isolated is a nSpots x 1 binary array giving 1 for
% well-isolated spots
%
% NB spots that can't be read in all rounds are discarded
% 
% This finds initial shifts between rounds using point cloud not by finding 
% the max correlation between images
%
% Optional second argument, AllBaseLocalYX, allows you to skip 
% first step if already have spots
%
% Optional third argument, allows you to skip PointCloudRegistrationStep.
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

%% variable naming conventions:
% spot subgroups:
% All: Any spot included in any tile (includes duplicates)
% nd: only spots whose anchor coordinate is in its home tile (no duplicates)
% Good: Spots for which could be read for all rounds

% coordinate frames or other info
% LocalYX: relative to home tile on the reference round
% LocalTile: number of home tile on the reference round
% GlobalYX: relative to the stitched image on the reference round
% RoundYX: relative to home tile after registration on each round
% RoundTile: number of home tile after registration on each round
% Isolated: binary number, saying if it is isolated
% SpotColors: the answer:

%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end
%% basic variables
rr = o.ReferenceRound;
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end

[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
o.TileCentre = 0.5*[o.TileSz+1,o.TileSz+1];
ImageRounds = setdiff(o.UseRounds,o.ReferenceRound);

if nargin>1 && sum(size(AllBaseLocalYX) ~= [nTiles,o.nBP,o.nRounds])
    error('AllBaseLocalYX has shape [%d, %d, %d] when it should be [%d, %d, %d]',...
        size(AllBaseLocalYX,1), size(AllBaseLocalYX,2), size(AllBaseLocalYX,3),...
        nTiles,o.nBP,o.nRounds);
end

if nargin>2 && (isempty(o.A) || isempty(o.D))
    error('Trying to skip PointCloudRegistration but o.A or o.D are empty');
end


%% get spot local coordinates in all colour channels and run PCR
if nargin<2
    AllBaseLocalYX = cell(nTiles,o.nBP, o.nRounds);    
    % loop through all tiles, finding PCR outputs
    fprintf('\nLocating spots in each colour channel of tile   ');
    for t=NonemptyTiles
        
        if t<10
            fprintf('\b%d',t);
        else
            fprintf('\b\b%d',t);
        end
        
        for r=ImageRounds
            % find spots whose home tile on round r is t
            % open file for this tile/round
            FileName = o.TileFiles{r,t};
            TifObj = Tiff(FileName);
            
            % now read in images for each base
            for b=o.UseChannels
                
                TifObj.setDirectory(o.FirstBaseChannel + b - 1);
                BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
                
                % find spots for base b on tile t - we will use this for point
                % cloud registration only, we don't use these detections to
                % detect colors, we read the colors off the
                % pointcloud-corrected positions of the spots detected in the reference round home tiles
                Spots = o.detect_spots(BaseIm,t,b,r);
                AllBaseLocalYX(t,b,r) = {Spots};
                
            end
            TifObj.close();
        end
    end
    fprintf('\n');
end

if nargin<3 || SkipRegistration == false
    %Add reference round to imaging spots
    if o.ReferenceRound~=o.AnchorRound
        AllBaseLocalYX(:,:,o.ReferenceRound) = o.RawLocalYX;
        o.PointCloudMethod = 2;     %Fit chromatic aberration separately if no anchor
        warning('o.ReferenceRound is different from o.AnchorRound so setting o.PointCloudMethod=2.');
    end
    %Correct o.RawLocalYX by removing spots that appear in multiple colour
    %channels
    o = o.remove_reference_duplicates(NonemptyTiles);
    
    %Save workspace at various stages so dont have to go back to the beginning
    %and often fails at PCR step.
    save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');
    
    %% Find initial shifts between rounds and then run PCR
    %Should have a initial search range for each round. If only provided one,
    %set all other rounds to the same range.
    %Do PCR before removing duplicates this time, because reference round
    %contains spots from multiple colour channels so need to adjust their
    %positions to take account of chromatic aberration.
    if size(o.FindSpotsSearch,1) == 1
        FindSpotsSearch = cell(o.nRounds,1);
        for r = o.UseRounds
            FindSpotsSearch{r} = o.FindSpotsSearch;
        end
        o.FindSpotsSearch = FindSpotsSearch;
        clear FindSpotsSearch
    end
    
    %Unless specified, set initial shift channel to be the one with largest
    %number of spots on tile/round with least spots.
    AllBaseSpotNo = cell2mat(cellfun(@size,AllBaseLocalYX,'uni',false));
    o.AllBaseSpotNo = AllBaseSpotNo(:,1:2:o.nBP*2,:);
    
    %If tile has too few spots, set to empty
    BadTiles = find(min(min(squeeze(sum(o.AllBaseSpotNo,2)),[],2),o.RawLocalNo)<o.OutlierMinScore);
    if ~isempty(BadTiles)
        warning('Setting tiles %d off as have too few spots',BadTiles);
        o.EmptyTiles(BadTiles) = 1;
        NonemptyTiles = find(~o.EmptyTiles)';
        if size(NonemptyTiles,2)==1
            NonemptyTiles = NonemptyTiles';
        end
    end
    
    %Allow for one bad tile so look at second worst tile
    a = o.AllBaseSpotNo(NonemptyTiles,:,ImageRounds);
    a(find(a==min(o.AllBaseSpotNo(NonemptyTiles,:,ImageRounds),[],1))) = inf;
    MinColorChannelSpotNo = min(min(a,[],1),[],3);
    clear a;
    if ~ismember(string(o.InitialShiftChannel),string(o.UseRounds))
        [~,o.InitialShiftChannel] = max(MinColorChannelSpotNo);
    end
    
    if MinColorChannelSpotNo(o.InitialShiftChannel)< o.MinSpots
        [BestSpots,BestChannel] = max(MinColorChannelSpotNo);
        if BestSpots >= o.MinSpots
            warning('Changing from Color Channel (%d) to Color Channel (%d) to find initial shifts.',o.InitialShiftChannel,BestChannel);
            o.InitialShiftChannel = BestChannel;
        else
            warning('Best Color Channel (%d) only has %d spots. Not enough for finding initial shifts. Will try to find some shifts with Fft method.'...
                ,BestChannel,BestSpots);
        end
    end
    
    o.D0 = zeros(nTiles,2,o.nRounds);
    Scores = zeros(nTiles,o.nRounds);
    ChangedSearch = zeros(o.nRounds,1);
    OutlierShifts = zeros(nTiles,2,o.nRounds);
    
    for t=NonemptyTiles
        for r = ImageRounds
            tic
            if min(o.AllBaseSpotNo(t,o.InitialShiftChannel,r),o.RawLocalNo(t))<o.OutlierMinScore
                warning('Tile %d, round %d, only has %d spots so using Fft method', t, r,...
                    min(o.AllBaseSpotNo(t,o.InitialShiftChannel,r),o.RawLocalNo(t)));
                [o.D0(t,:,r), Scores(t,r)] = o.get_Fft_shift_single(t,r,o.InitialShiftChannel,...
                    t,o.ReferenceRound,o.ReferenceChannel,'FindSpots');
                tChangedSearch = 0;
            else
                [o.D0(t,:,r), Scores(t,r),tChangedSearch] = o.get_initial_shift2(AllBaseLocalYX{t,o.InitialShiftChannel,r},...
                    vertcat(o.RawLocalYX{t,:}), o.FindSpotsSearch{r},'FindSpots');
            end
            toc
            ChangedSearch(r) = ChangedSearch(r)+tChangedSearch;
            
            fprintf('Tile %d, shift from reference round (%d) to round %d: [%d %d], score %d \n',...
                t, o.ReferenceRound,r, o.D0(t,:,r),Scores(t,r));
            
            %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
            if t == 3 || (mod(ChangedSearch(r),2) == 0) && (ChangedSearch(r)>0)
                o = o.GetNewSearchRange_FindSpots(t,r);
            end
        end
    end
    
    for r = ImageRounds
        [o.D0(NonemptyTiles,:,r), OutlierShifts(NonemptyTiles,:,r)] = ...
            o.AmendShifts(o.D0(NonemptyTiles,:,r),Scores(NonemptyTiles,r),'FindSpots');
    end
    
    o.FindSpotsInfo.Scores = Scores;
    o.FindSpotsInfo.ChangedSearch = ChangedSearch;
    o.FindSpotsInfo.Outlier = OutlierShifts;
    
    save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');
    
    if o.PointCloudMethod==1
        [o,x] = o.PointCloudRegister6(AllBaseLocalYX, o.RawLocalYX, 1, nTiles);
    elseif o.PointCloudMethod==2
        [o,x] = o.PointCloudRegister2(AllBaseLocalYX, o.RawLocalYX, 1, nTiles);
    else
        error('Require o.PointCloudMethod to be 1 or 2');
    end
    %Reference round coordinates are adjusted as to take account of chromatic
    %aberration.
    o.RawLocalYX = x;
    %o = AmendPCR(o);
    o = PCR_analysis(o,AllBaseLocalYX,o.RawLocalYX);
    
    clearvars x;
    
    save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');
    
    if strcmpi(o.PcImageMatchesThresh, 'auto')
        o.PcImageMatchesThresh = length(NonemptyTiles)*7;
    end
    if ~isnumeric(o.MinPCMatchFract) || o.MinPCMatchFract>=1
        o.MinPCMatchFract = 0.1;
    end
    if ~isnumeric(o.MaxPCMatchThresh)
        o.MaxPCMatchThresh = 500;
    end
    
    nMatches = o.nMatches(NonemptyTiles,:,:);
    AllBaseSpotNo = o.AllBaseSpotNo(NonemptyTiles,:,:);
    PcMatchThresh = o.MinPCMatchFract*AllBaseSpotNo;
    PcMatchThresh(PcMatchThresh>o.MaxPCMatchThresh) = o.MaxPCMatchThresh;
    nBadRegImages = length(nMatches(nMatches<PcMatchThresh));
    if nBadRegImages>o.PcImageMatchesThresh
        ErrorFile = fullfile(o.OutputDirectory, 'oFindSpots-Error_with_PointCloudRegistration');
        save(ErrorFile, 'o', '-v7.3');
        error(['%d/%d images have nMatches<o.MinPCMatchFract*o.AllBaseSpotNo, where o.MinPCMatchFract =  %.2f'...
            '\nThis exceeds threshold of o.PcImageMatchesThresh = %d.'...
            '\nProgress up to this point saved as:\n%s.mat'],...
            nBadRegImages,nTiles*o.nBP*o.nRounds,o.MinPCMatchFract,o.PcImageMatchesThresh,ErrorFile);
    end
end

%% now make array of global coordinates
nAll = sum(sum(cellfun(@numel, o.RawIsolated)));

AllGlobalYX = zeros(nAll,2);
AllLocalYX = zeros(nAll,2);
AllIsolated = zeros(nAll,1);
OriginalTile = zeros(nAll,1);
AllOriginalChannel = zeros(nAll,1);     %Keep track of which channel each spot from

ind = 1;
for t=NonemptyTiles
    ChannelSizes = cellfun(@numel, o.RawIsolated(t,:));
    nMySpots = sum(ChannelSizes);
    AllGlobalYX(ind:ind+nMySpots-1,:) = vertcat(o.RawLocalYX{t,:})+o.TileOrigin(t,:,rr);
    AllLocalYX(ind:ind+nMySpots-1,:) = vertcat(o.RawLocalYX{t,:});
    AllIsolated(ind:ind+nMySpots-1) = vertcat(o.RawIsolated{t,:});
    OriginalTile(ind:ind+nMySpots-1) = t;
    AllOriginalChannel(ind:ind+nMySpots-1) = repelem(o.UseChannels,ChannelSizes(o.UseChannels))';
    ind = ind+nMySpots;
end

if o.Graphics
    figure(1001)
    colors = colormap(lines(o.nBP));
    hold on
    for b=o.ReferenceSpotChannels
        ToPlot = AllOriginalChannel == b;
        plot(AllGlobalYX(ToPlot,2), AllGlobalYX(ToPlot,1), '.', 'markersize', 4,'Color',colors(b,:));        
    end
    title('All global coords including duplicates');
    leg = legend(o.bpLabels(o.ReferenceSpotChannels),'Location','northwest');
    title(leg,'Color Channel');
    hold off
    %set(gca, 'YDir', 'reverse');
end

%% now remove duplicates by keeping only spots detected on their home tile

[AllLocalTile, ~] = which_tile(AllGlobalYX, o.TileOrigin(:,:,rr), o.TileSz);
NotDuplicate = (AllLocalTile==OriginalTile);
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndLocalYX = AllLocalYX(NotDuplicate,:);
ndIsolated = AllIsolated(NotDuplicate,:);
ndLocalTile = AllLocalTile(NotDuplicate,:);
ndOriginalChannel = AllOriginalChannel(NotDuplicate,:);

nnd = sum(NotDuplicate);

if o.Graphics
    figure(1002); clf
    colors = colormap(lines(o.nBP));
    hold on
    for b=o.ReferenceSpotChannels
        ToPlot = ndOriginalChannel == b;
        plot(ndGlobalYX(ToPlot,2), ndGlobalYX(ToPlot,1), '.', 'markersize', 4,'Color',colors(b,:));        
    end
    title('Global coords without duplicates');
    leg = legend(o.bpLabels(o.ReferenceSpotChannels),'Location','northwest');
    title(leg,'Color Channel');
    hold off
    drawnow;
    %set(gca, 'YDir', 'reverse');
end


%% decide which tile to read each spot off in each round. 
% They are read of home tile if possible (always possible in ref round)
% in other rounds might have to be a NWSE neighbor - but never a diagonal
% neighbor
% ndRoundTile(s,r) stores appropriate tile for spot s on round r
% ndRoundYX(s,:,r) stores YX coord on this tile

%Compute approx new shifts from D matrices
YXRoundTileShifts = permute(squeeze(mean(o.D(3,:,:,1:o.nRounds,:),5)),[2,1,3]);
o.TileOrigin(:,:,1:o.nRounds) =  o.TileOrigin(:,:,rr) - YXRoundTileShifts;  

[ndRoundTile,ndRoundYX] = o.get_SpotTileEachRound(ndGlobalYX,ndLocalTile);

%% loop through all tiles, finding spot colors
[ndSpotColors,ndPointCorrectedLocalYX] = get_spot_colors(o,ndLocalYX,ndLocalTile,...
    ndRoundTile,AllBaseLocalYX,o.nMatches,o.Error);

%% now find those that were detected in all tiles
ndSpotColorsToUse = ndSpotColors(:,o.UseChannels,o.UseRounds);
Good = all(isfinite(ndSpotColorsToUse(:,:)),2);
GoodGlobalYX = ndGlobalYX(Good,:);
GoodSpotColors = ndSpotColors(Good,:,:);
GoodLocalTile = ndLocalTile(Good);
GoodIsolated = ndIsolated(Good);

save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX',...
    'Good', 'ndGlobalYX', 'ndSpotColors', 'ndLocalTile','ndIsolated',...
    'ndPointCorrectedLocalYX','ndRoundYX','ndRoundTile');

%% plot those that were found and those that weren't
if o.Graphics
    figure(10032); clf; hold on; set(gca, 'color', 'k');
    plot(ndGlobalYX(Good,2), ndGlobalYX(Good,1), 'b.', 'markersize', 1);
    plot(ndGlobalYX(~Good,2), ndGlobalYX(~Good,1), 'r.', 'markersize', 1);
    legend({'resolved', 'unresolved'}, 'color', [.6 .6 .6]);
    % now put on edges
    SquareX1 = [0, 0, o.TileSz];
    SquareY1 = [o.TileSz, 0, 0];
    SquareX2 = [o.TileSz, o.TileSz, 0];
    SquareY2 = [0, o.TileSz, o.TileSz];

    SquareColors = hsv2rgb([(1:o.nRounds)'/o.nRounds, [.5, .6] .*ones(o.nRounds,1)]);
    SquareColors(o.ReferenceRound,:)=1.0;
    for r=o.UseRounds
        for t=NonemptyTiles
            MyOrigin = o.TileOrigin(t,:,r);
            plot(SquareX1 + MyOrigin(2), SquareY1 + MyOrigin(1),...
                '--', 'Color', SquareColors(r,:));
            plot(SquareX2 + MyOrigin(2), SquareY2 + MyOrigin(1),...
                ':', 'Color', SquareColors(r,:));

            text(MyOrigin(2), MyOrigin(1),...
                sprintf('T%d r%d', t, r), 'color', SquareColors(r,:)); 
        end
    end
    
    
    %set(gca, 'YDir', 'reverse');
end
       

%%
o.dpSpotGlobalYX = GoodGlobalYX;
o.dpSpotColors = GoodSpotColors;          
o.dpSpotIsolated = logical(GoodIsolated);
o.dpLocalTile = int16(GoodLocalTile);
