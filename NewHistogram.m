Dist = 25;      %Get background from within Dist pixels of anchor spots
nCounts = 15;   %Only add pixel to background pixels if is within Dist of nCounts anchor spots
rr = o.ReferenceRound;
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';

%% Find all coordinates within Dist of anchor spots
AllLocalYX = zeros(o.TileSz^2,2);
AllLocalYX(:,1) = repelem(1:o.TileSz,1,o.TileSz);
AllLocalYX(:,2) = repmat(1:o.TileSz,1,o.TileSz);


BackgroundLocalYX = cell(nTiles,1);

for t = NonemptyTiles
    Idx = rangesearch(AllLocalYX,o.RawLocalYX{t,o.ReferenceChannel},Dist);
    [UniqueIdx,~,iUniqueIdx] = unique(horzcat(Idx{:})');
    UniqueIdx_counts = accumarray(iUniqueIdx,1);
    BackgroundLocalYX{t} = AllLocalYX(UniqueIdx(UniqueIdx_counts>nCounts),:);
end


%% now make array of global coordinates
nAll = sum(sum(cellfun(@numel, BackgroundLocalYX)/2));

AllGlobalYX = zeros(nAll,2);
AllLocalYX = zeros(nAll,2);
AllIsolated = zeros(nAll,1);
OriginalTile = zeros(nAll,1);
AllOriginalChannel = zeros(nAll,1);     %Keep track of which channel each spot from

ind = 1;
for t=NonemptyTiles
    nMySpots = length(BackgroundLocalYX{t});
    AllGlobalYX(ind:ind+nMySpots-1,:) = vertcat(BackgroundLocalYX{t})+o.TileOrigin(t,:,rr);
    AllLocalYX(ind:ind+nMySpots-1,:) = vertcat(BackgroundLocalYX{t});
    OriginalTile(ind:ind+nMySpots-1) = t;
    ind = ind+nMySpots;
end
if o.Graphics
    figure(1001); clf
    plot(AllGlobalYX(:,2), AllGlobalYX(:,1), '.', 'markersize', 4);        
    title('All global coords including duplicates');
    drawnow;
    %set(gca, 'YDir', 'reverse');
end

%% now remove duplicates by keeping only spots detected on their home tile

[AllLocalTile, ~] = which_tile(AllGlobalYX, o.TileOrigin(:,:,rr), o.TileSz);
NotDuplicate = (AllLocalTile==OriginalTile);
ndGlobalYX = AllGlobalYX(NotDuplicate,:);
ndLocalYX = AllLocalYX(NotDuplicate,:);
ndLocalTile = AllLocalTile(NotDuplicate,:);
nnd = sum(NotDuplicate);

if o.Graphics
    figure(1002); clf
    plot(ndGlobalYX(:,2), ndGlobalYX(:,1), '.', 'markersize', 4);        
    title('Global coords without duplicates');
    drawnow;
    %set(gca, 'YDir', 'reverse');
end

%% decide which tile to read each spot off in each round. 
% They are read of home tile if possible (always possible in ref round)
% in other rounds might have to be a NWSE neighbor - but never a diagonal
% neighbor
% ndRoundTile(s,r) stores appropriate tile for spot s on round r
% ndRoundYX(s,:,r) stores YX coord on this tile

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

%% Read off intensities from all imaging rounds
ndLocalYX = [ndLocalYX-o.TileCentre,ones(nnd,1)];
ndSpotColors = nan(nnd, o.nBP, o.nRounds);
load([o.OutputDirectory,'/FindSpotsWorkspace.mat'],'AllBaseLocalYX');
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
                    MyPointCorrectedYX = MyLocalYX*o.D(:,:,t,r,b)+o.TileCentre;
                    MyPointCorrectedYX = round(MyPointCorrectedYX);
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                else
                    [MyPointCorrectedYX, Error, nMatches] = o.different_tile_transform(AllBaseLocalYX,o.RawLocalYX, ...
                        MyLocalYX,t,t2,r,b);
                    fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                        t, t2, r, b,  nMatches, o.RawLocalNo(t2), Error);
                    if nMatches<o.MinPCMatchFract*o.AllBaseSpotNo(t,b,r) || isempty(nMatches)
                        continue;
                    end
                    ndSpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                end
               
            end    
        end
        TifObj.close();       
    end
end
fprintf('\n');
ndSpotColorsToUse = ndSpotColors(:,o.UseChannels,o.UseRounds);
Good = all(isfinite(ndSpotColorsToUse(:,:)),2);
GoodSpotColors = ndSpotColors(Good,:,:);


%% Make histograms
NewHistCounts = zeros(length(o.HistValues),o.nBP,o.nRounds);
for b=1:o.nBP
    for r=1:o.nRounds
        NewHistCounts(:,b,r) = histc(GoodSpotColors(:,b,r),o.HistValues);
    end
end

%Plot histograms to make sure they are smooth
%Avoid ExtraRounds as only need histograms for the 7 rounds used to
%define genes

if o.Graphics
    figure(43291);
    index = 1;
    for r=1:o.nRounds
        for b=1:o.nBP
            nPixels = sum(o.HistCounts(:,b,r));
            nPixelsNew = sum(NewHistCounts(:,b,r));
            subplot(o.nRounds,o.nBP,index)
            histogram('BinEdges',[o.HistValues-0.5,max(o.HistValues)+0.5],'BinCounts',o.HistCounts(:,b,r)/nPixels,'DisplayStyle','stairs');
            %xlim([-1000,1000]);
            ylim([0,max(o.HistCounts(:,b,r))/nPixels]);
            hold on
            histogram('BinEdges',[o.HistValues-0.5,max(o.HistValues)+0.5],'BinCounts',NewHistCounts(:,b,r)/nPixelsNew,'DisplayStyle','stairs');
            hold off
            if b==4
                title(strcat('Round ',num2str(r)));
            end
            index = index+1;
        end
    end
end