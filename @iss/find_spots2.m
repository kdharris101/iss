function o = find_spots2(o)        
% o = o.find_spots2;
%
% finds spots in all tiles using the reference channel, removes
% duplicates in overlap regions and returns nSpots x 2 array o.SpotGlobalYX of
% coordinates in global frame
% 
% Looks up colors from apporpriate files and makes nSpots x nBP x nRounds
% array o.SpotColors
%
% o.Isolated is a nSpots x 1 binary array giving 1 for
% well-isolated spots
%
% NB spots that can't be read in all rounds are discarded
% 
% This finds initial shifts between rounds using point cloud not by finding 
% the max correlation between images
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

%% basic variables
rr = o.ReferenceRound;
Tiles = find(~o.EmptyTiles)';

[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;

    
%% now make array of global coordinates
AllIsolated = logical(vertcat(o.RawIsolated{:})); % I HATE MATLAB - for converting logical to doubles for no reason
nAll = length(AllIsolated);

AllGlobalYX = zeros(nAll,2);
AllLocalYX = zeros(nAll,2);

ind = 1;
for t=Tiles
    MySpots = o.RawLocalYX{t};
    nMySpots = size(MySpots, 1);
    AllGlobalYX(ind:ind+nMySpots-1,:) = bsxfun(@plus, MySpots, o.TileOrigin(t,:,rr));
    AllLocalYX(ind:ind+nMySpots-1,:) = MySpots;
    OriginalTile(ind:ind+nMySpots-1) = t;
    ind = ind+nMySpots;
end
if o.Graphics
    figure(1001)
    plot(AllGlobalYX(:,2), AllGlobalYX(:,1), '.', 'markersize', 1);
    title('All global coords including duplicates');
    %set(gca, 'YDir', 'reverse');
end

%% now remove duplicates by keeping only spots detected on their home tile

[AllLocalTile, ~] = which_tile(AllGlobalYX, o.TileOrigin(:,:,rr), o.TileSz);
NotDuplicate = (AllLocalTile==OriginalTile');
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndLocalYX = AllLocalYX(NotDuplicate,:);
ndIsolated = AllIsolated(NotDuplicate,:);
ndLocalTile = AllLocalTile(NotDuplicate,:);

nnd = sum(NotDuplicate);

if o.Graphics
    figure(1002); clf
    plot(ndGlobalYX(:,2), ndGlobalYX(:,1), '.', 'markersize', 1);
    title('Global coords without duplicates');
    drawnow;
    %set(gca, 'YDir', 'reverse');
end




%% get spot local coordinates in all colour channels and run PCR

%Specify which rounds/colour channels to use (default is all)
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end

if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

AllBaseLocalYX = cell(nTiles,o.nBP, o.nRounds);

SE = fspecial('disk', o.SmoothSize);

% loop through all tiles, finding PCR outputs
fprintf('\nLocating spots in each colour channel of tile   ');
for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    
    if t<10
        fprintf('\b%d',t);
    else
        fprintf('\b\b%d',t);
    end 
    
    [y, x] = ind2sub([nY nX], t);

    for r=o.UseRounds
        % find spots whose home tile on round r is t      
        % open file for this tile/round
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
                
        % now read in images for each base
        for b=o.UseChannels             

            TifObj.setDirectory(o.FirstBaseChannel + b - 1);            
            BaseIm = TifObj.read();
            
            % find spots for base b on tile t - we will use this for point
            % cloud registration only, we don't use these detections to
            % detect colors, we read the colors off the
            % pointcloud-corrected positions of the spots detected in the reference round home tiles  
            CenteredSpots = o.detect_spots(BaseIm) - [o.TileSz/2,o.TileSz/2];
            AllBaseLocalYX(t,b,r) = {CenteredSpots};

        end
        TifObj.close();
    end      
end
fprintf('\n');

%Save workspace at various stages so dont have to go back to the beginning
%and often fails at PCR step.
save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');

%Should have a initial search range for each round. If only provided one,
%set all other rounds to the same range.
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
o.AllBaseSpotNo = AllBaseSpotNo(:,1:2:o.nRounds*2,:);
MinColorChannelSpotNo = min(min(o.AllBaseSpotNo),[],3);
if ~ismember(string(o.InitialShiftChannel),string(o.UseRounds))
    [~,o.InitialShiftChannel] = max(MinColorChannelSpotNo);
end

if MinColorChannelSpotNo(o.InitialShiftChannel)< o.MinSpots
    [BestSpots,BestChannel] = max(MinColorChannelSpotNo);
    if BestSpots >= o.MinSpots
        warning('Changing from Color Channel (%d) to Color Channel (%d) to find initial shifts.',o.InitialShiftChannel,BestChannel);
        o.InitialShiftChannel = BestChannel;
    else
        error('Best Color Channel (%d) only has %d spots. Not enough for finding initial shifts. Consider reducing o.DetectionThresh.'...
            ,BestChannel,BestSpots);
    end
end

o.D0 = zeros(nTiles,2,o.nRounds);
Scores = zeros(nTiles,o.nRounds);
ChangedSearch = zeros(o.nRounds,1);
OutlierShifts = zeros(nTiles,2,o.nRounds);

for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    for r = o.UseRounds
        tic
        [o.D0(t,:,r), Scores(t,r),tChangedSearch] = o.get_initial_shift2(AllBaseLocalYX{t,o.InitialShiftChannel,r},...
            o.RawLocalYX{t}, o.FindSpotsSearch{r},'FindSpots');
        toc
        ChangedSearch(r) = ChangedSearch(r)+tChangedSearch;
        
        fprintf('Tile %d, shift from anchor round to round %d: [%d %d], score %d \n', t, r, o.D0(t,:,r),...
            Scores(t,r));
        
        %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
        if t == 3 || (mod(ChangedSearch(r),2) == 0) && (ChangedSearch(r)>0)
            o = o.GetNewSearchRange_FindSpots(t,r);
        end
    end
end


for r = o.UseRounds
    [o.D0(:,:,r), OutlierShifts(:,:,r)] = o.AmendShifts(o.D0(:,:,r),Scores(:,r),'FindSpots');
end

o.FindSpotsInfo.Scores = Scores;
o.FindSpotsInfo.ChangedSearch = ChangedSearch;
o.FindSpotsInfo.Outlier = OutlierShifts;

save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');

o = o.PointCloudRegister(AllBaseLocalYX, o.RawLocalYX, eye(2), nTiles);

save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');

%% decide which tile to read each spot off in each round. 
% They are read of home tile if possible (always possible in ref round)
% in other rounds might have to be a NWSE neighbor - but never a diagonal
% neighbor
% ndRoundTile(s,r) stores appropriate tile for spot s on round r
% ndRoundYX(s,:,r) stores YX coord on this tile

o.TileOrigin(:,:,1:o.nRounds) =  o.TileOrigin(:,:,rr) - o.D;

ndRoundTile = nan(nnd,o.nRounds);
ndRoundYX = nan(nnd,2,o.nRounds);

PossNeighbs = [-1 -nY 1 nY 0]; % NWSE then same tile - same will have priority by being last

for r=o.UseRounds
    fprintf('Finding appropriate tiles for round %d\n', r);
    
    for n = PossNeighbs
        % find origins of each tile's neighbor, NaN if not there
        NeighbTile = (1:nTiles)+n;
        NeighbOK = (NeighbTile>=1 & NeighbTile<=nTiles);
        NeighbOrigins = nan(nTiles,2);
        NeighbOrigins(NeighbOK,:) = round(o.TileOrigin(NeighbTile(NeighbOK),:,r));
        
        % now for each spot see if it is inside neighbor's tile area
        SpotsNeighbOrigin = NeighbOrigins(ndLocalTile,:);
        SpotsInNeighbTile = all(ndGlobalYX>=SpotsNeighbOrigin+1+o.ExpectedAberration...
            & ndGlobalYX<=SpotsNeighbOrigin+o.TileSz-o.ExpectedAberration, 2);
        
        % for those that were in set this to be its neighbor
        ndRoundTile(SpotsInNeighbTile,r) = NeighbTile(ndLocalTile(SpotsInNeighbTile));    
    end
    
    % compute YX coord
    HasTile = isfinite(ndRoundTile(:,r));
    ndRoundYX(HasTile,:,r) = ndGlobalYX(HasTile,:) - round(o.TileOrigin(ndRoundTile(HasTile,r),:,r));
    
end

%% loop through all tiles, finding spot colors

ndSpotColors = nan(nnd, o.nBP, o.nRounds);
ndPointCorrectedLocalYX = nan(nnd, 2, o.nRounds, o.nBP);

for t=1:nTiles
    if o.EmptyTiles(t); continue; end
    [y, x] = ind2sub([nY nX], t);
   
    for r=o.UseRounds         
        % find spots whose home tile on round r is t
        MySpots = (ndRoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        
        % open file for this tile/round
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = ndLocalTile(ndRoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        fprintf('\nRef round home tiles for spots in t%d at (%2d, %2d), r%d: ', t, y, x, r);
        for i=MyRefTiles(:)'
            fprintf('t%d, %d spots; ', i, sum(RefRoundHomeTiles==i));
        end
        fprintf('\n');        
        
        
        % now read in images for each base
        for b=o.UseChannels               %No 0 as trying without using anchor

            
            TifObj.setDirectory(o.FirstBaseChannel + b - 1);
            BaseIm = TifObj.read();
            
            if o.SmoothSize
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
            else
                BaseImSm = BaseIm;
            end
            
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (ndRoundTile(:,r)==t & ndLocalTile==t2);
                CenteredMyLocalYX = ndLocalYX(MyBaseSpots,:) - [o.TileSz/2,o.TileSz/2];
                
                if t == t2
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  o.nMatches(t,b,r), size(o.RawLocalYX{t2},1), o.Error(t,b,r));
                    if o.nMatches(t,b,r)<o.MinPCMatches || isempty(o.nMatches(t,b,r))
                        continue;
                    end
                    CenteredMyPointCorrectedYX = (o.A(:,:,b)*(CenteredMyLocalYX + o.D(t,:,r))')';
                    MyPointCorrectedYX = round(CenteredMyPointCorrectedYX + [o.TileSz/2,o.TileSz/2]);
                    ndPointCorrectedLocalYX(MyBaseSpots,:,r,b) = MyPointCorrectedYX;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                else
                    [MyPointCorrectedYX, Error, nMatches] = o.different_tile_transform(AllBaseLocalYX,o.RawLocalYX, ...
                        CenteredMyLocalYX,t,t2,r,b);
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  nMatches, size(o.RawLocalYX{t2},1), Error);
                    if nMatches<o.MinPCMatches || isempty(nMatches)
                        continue;
                    end
                    ndPointCorrectedLocalYX(MyBaseSpots,:,r,b) = MyPointCorrectedYX;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                end
               
            end    
        end
        TifObj.close();       
    end
end
fprintf('\n');

%% now find those that were detected in all tiles
ndSpotColorsToUse = ndSpotColors(:,o.UseChannels,o.UseRounds);
Good = all(isfinite(ndSpotColorsToUse(:,:)),2);
GoodGlobalYX = ndGlobalYX(Good,:);
GoodSpotColors = ndSpotColors(Good,:,:);
GoodLocalTile = ndLocalTile(Good);
GoodIsolated = ndIsolated(Good);

save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX',...
    'Good', 'ndGlobalYX', 'ndSpotColors', 'ndLocalTile','ndIsolated','ndPointCorrectedLocalYX','ndRoundYX','ndRoundTile');

%% plot those that were found and those that weren't
if o.Graphics
    PlotScale = 1;
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
        for t=Tiles
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
       


%% sanity check
plsz = 7;
if o.Graphics ==2
    
    GoodRoundYX = ndRoundYX(Good,:,:);
    GoodRoundTile = ndRoundTile(Good,:);
    GoodCorrectedYX = ndPointCorrectedLocalYX(Good,:,:,:);

    roi = o.FindSpotsRoi;
    PlotSpots = find(GoodGlobalYX(:,1)>roi(1) & GoodGlobalYX(:,1)<roi(2) & GoodGlobalYX(:,2)>roi(3) & GoodGlobalYX(:,2)<roi(4));
    
    for s=(PlotSpots(:))' %PlotSpots(randperm(length(PlotSpots)))'
        figure(91); clf
        for r=o.UseRounds
            t=GoodRoundTile(s,r);

            fprintf('Spot %d, round %d, tile %d: y=%d, x=%d\n', s, r, t, GoodRoundYX(s,1,r), GoodRoundYX(s,2,r));

            Ylegends = {'Anchor', o.bpLabels{:}};
            for b=o.UseChannels
                
                      
%                 if b==0                    
%                     y0 = GoodRoundYX(s,1,r);
%                     x0 = GoodRoundYX(s,2,r);
%                 else
%                     y0 = GoodCorrectedYX(s,1,r,b);
%                     x0 = GoodCorrectedYX(s,2,r,b);
%                 end
                y0 = GoodCorrectedYX(s,1,r,b);
                x0 = GoodCorrectedYX(s,2,r,b);
                if ~isfinite(x0) || ~isfinite(y0)
                    continue;
                end
                y1 = max(1,y0 - plsz);
                y2 = min(o.TileSz,y0 + plsz);
                x1 = max(1,x0 - plsz);
                x2 = min(o.TileSz,x0 + plsz);
           
                
                BaseIm = imread(o.TileFiles{r,t}, b, 'PixelRegion', {[y1 y2], [x1 x2]});
                if o.SmoothSize
                    BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
                else
                    BaseImSm = BaseIm;
                end

                subplot(o.nBP+1, o.nRounds, (b)*o.nRounds + r)
                imagesc([x1 x2], [y1 y2], BaseImSm); hold on
                axis([x0-plsz, x0+plsz, y0-plsz, y0+plsz]);
                plot(xlim, [y0 y0], 'w'); plot([x0 x0], ylim, 'w');
                caxis([0 o.DetectionThresh*2]);
                if r==1; ylabel(Ylegends{b+1}); end
                colorbar;
                
                title(sprintf('Round %d, Base %d, Tile %d', r, b, t));
                drawnow
            end
        end
        fprintf('\n');
        figure(92); clf
        imagesc(sq(GoodSpotColors(s,:,:)));
        set(gca, 'ytick', 1:5); set(gca, 'yticklabel', {'Anchor', o.bpLabels{:}});
        %caxis([0 o.DetectionThresh*2]);
%         fprintf('local YX = (%f, %f) screen YX = (%f, %f) Called as %s, %s, quality %f\n', ...
%             GoodRoundYX(s,1), GoodRoundYX(s,2), GoodGlobalYX(s,1)/4, GoodGlobalYX(s,2)/4, ...
%             GoodCodes{s}, GoodGenes{s}, GoodMaxScore(s));
        figure(1003); hold on
        squarex = [-1 1 1 -1 -1]*plsz; squarey = [-1 -1 1 1 -1]*plsz;
        h = plot(GoodGlobalYX(s,2)+squarex, GoodGlobalYX(s,1)+squarey, 'g');
        pause;
        delete(h);
    end
end



%%
o.SpotGlobalYX = GoodGlobalYX;
o.cSpotColors = GoodSpotColors;          
%o.cAnchorIntensities = squeeze(GoodSpotColors(:,1,:));
o.cSpotIsolated = GoodIsolated;
