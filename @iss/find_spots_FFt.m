function o = find_spots_FFt(o)        
% o = o.find_spots_FFt;
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


%% get spot local coordinates in all colour channels and find initial shifts

AllBaseLocalYX = cell(nTiles,o.nBP, o.nRounds);
ImageRounds = setdiff(o.UseRounds,o.ReferenceRound);
o.D0 = zeros(nTiles,2,o.nRounds);
o.cc = zeros(nTiles,o.nRounds);
% loop through all tiles, finding PCR outputs
for t=NonemptyTiles
    
    [y, x] = ind2sub([nY nX], t);
    
    %Reload anchor image to find initial shift
    AnchorIm = int32(imread(o.TileFiles{rr,y,x}, o.AnchorChannel)-o.TilePixelValueShift);
    if o.SmoothSize
        AnchorImSm = uint16(imfilter(AnchorIm, SE));
    else
        AnchorImSm = uint16(AnchorIm);
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
            if b == o.InitialShiftChannel
                %For chosen channel, find initial shift
                %Deliberately set negative pixels to 0, easier to register
                [o.D0(t,:,r), o.cc(t,r)] = o.ImRegFft2_FindSpots(uint16(BaseIm),AnchorImSm, 0, o.RegMinSize);
                fprintf('Tile %d, shift from reference round (%d) to round %d: [%d %d], score %d \n',...
                        t, o.ReferenceRound, r, o.D0(t,:,r),o.cc(t,r));
            end
            
        end
        TifObj.close();
    end      
end
fprintf('\n');

OutlierShifts = zeros(nTiles,2,o.nRounds);
for r = ImageRounds
    [o.D0(NonemptyTiles,:,r), OutlierShifts(NonemptyTiles,:,r)] = ...
        o.AmendShifts(o.D0(NonemptyTiles,:,r),o.OutlierMinScore-o.cc(NonemptyTiles,r),'FindSpots');
end
o.FindSpotsInfo.Scores = o.cc;
o.FindSpotsInfo.Outlier = OutlierShifts;

%Add reference round to imaging spots
if o.ReferenceRound~=o.AnchorRound
    AllBaseLocalYX(:,:,o.ReferenceRound) = o.RawLocalYX;
end
%Correct o.RawLocalYX by removing spots that appear in multiple colour
%channels
o = o.remove_reference_duplicates(NonemptyTiles);

%Save workspace at various stages so dont have to go back to the beginning
%and often fails at PCR step.
save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');

%% Run PCR

%Unless specified, set initial shift channel to be the one with largest
%number of spots on tile/round with least spots.
AllBaseSpotNo = cell2mat(cellfun(@size,AllBaseLocalYX,'uni',false));
o.AllBaseSpotNo = AllBaseSpotNo(:,1:2:o.nRounds*2,:);

[o,x] = o.PointCloudRegister2(AllBaseLocalYX, o.RawLocalYX, 1, nTiles);
%Reference round coordinates are adjusted as to take account of chromatic
%aberration.
o.RawLocalYX = x;
clearvars x;

save(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'o', 'AllBaseLocalYX');

if strcmpi(o.PcImageMatchesThresh, 'auto')
    o.PcImageMatchesThresh = 1.5*length(NonemptyTiles);
end

if ~isnumeric(o.MinPCMatchFract) || o.MinPCMatchFract>=1
    o.MinPCMatchFract = 0.1;
end

nBadRegImages = length(o.nMatches(o.nMatches<o.MinPCMatchFract*o.AllBaseSpotNo));
if nBadRegImages>o.PcImageMatchesThresh
    ErrorFile = fullfile(o.OutputDirectory, 'oFindSpots-Error_with_PointCloudRegistration');
    save(ErrorFile, 'o', '-v7.3');
    error(['%d/%d images have nMatches<o.MinPCMatchFract*o.AllBaseSpotNo, where o.MinPCMatchFract =  %.2f'...
        '\nThis exceeds threshold of o.PcImageMatchesThresh = %d.'...
        '\nProgress up to this point saved as:\n%s.mat'],...
        nBadRegImages,nTiles*o.nBP*o.nRounds,o.MinPCMatchFract,o.PcImageMatchesThresh,ErrorFile);
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
YXRoundTileShifts = permute(squeeze(o.D(3,:,:,1:o.nRounds)),[2,1,3]);
%Below used to be TileOrigin-YXShifts, but I think it should be plus, as
%when transfer from reference round to imaging round, you add shift. Hence
%when you transfer origin from reference round to imaging round, you should
%also add shift. Not sure though. !!!!!!!POSSIBLY BUG HERE!!!!!!!
o.TileOrigin(:,:,1:o.nRounds) =  o.TileOrigin(:,:,rr) + YXRoundTileShifts;  

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
ndLocalYX = [ndLocalYX-o.TileCentre,ones(nnd,1)];
ndSpotColors = nan(nnd, o.nBP, o.nRounds);
ndPointCorrectedLocalYX = nan(nnd, 2, o.nRounds, o.nBP);

for t=NonemptyTiles
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
            BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
            
            if o.SmoothSize
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
            else
                BaseImSm = BaseIm;
            end
            
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (ndRoundTile(:,r)==t & ndLocalTile==t2);
                MyLocalYX = ndLocalYX(MyBaseSpots,:);
                
                if t == t2
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  o.nMatches(t,b,r), o.RawLocalNo(t2), o.Error(t,b,r));
                    if o.nMatches(t,b,r)<o.MinPCMatchFract*o.AllBaseSpotNo(t,b,r) || isempty(o.nMatches(t,b,r))
                        warning('Tile %d, channel %d, round %d has %d point cloud matches, which is below the threshold of %d.',...
                            t,b,r,o.nMatches(t,b,r),o.MinPCMatchFract*o.AllBaseSpotNo(t,b,r));
                    end
                    MyPointCorrectedYX = o.A(b)*(MyLocalYX*o.D(:,:,t,r))+o.TileCentre;
                    MyPointCorrectedYX = round(MyPointCorrectedYX);
                    ndPointCorrectedLocalYX(MyBaseSpots,:,r,b) = MyPointCorrectedYX;
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                else
                    [MyPointCorrectedYX, Error, nMatches] = o.different_tile_transform(AllBaseLocalYX,o.RawLocalYX, ...
                        MyLocalYX,t,t2,r,b);
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  nMatches, o.RawLocalNo(t2), Error);
                    if nMatches<o.MinPCMatchFract*o.AllBaseSpotNo(t,b,r) || isempty(nMatches)
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
o.SpotGlobalYX = GoodGlobalYX;
o.cSpotColors = GoodSpotColors;          
%o.cAnchorIntensities = squeeze(GoodSpotColors(:,1,:));
o.cSpotIsolated = logical(GoodIsolated);
