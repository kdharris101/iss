function o=register2(o)
% SEEMS SLOWER THAN register.m SO WOULDN'T RECOMMEND USING
% o=iss_register2(o)
% 
% register images based on tile files
% creates arrays o.RefPos(y, x): origin of tile (y,x) in pixels on
% reference round relative to global coordinate frame
% 
% o.RelativePos(r, 1:2, t1, t2): origin of tile t2 on 
% reference round minus origin of tile t1 round r. In other words,
% Im2(x;rr) = Im1(x + RelativePos; r). Nan if not defined
% t1 and t2 are linear indices for the tile (y,x)
%
% also makes a global DAPI image 
% 
% This finds shifts between tiles using point cloud not by finding the max
% correlation between images. IN 2D, THIS APPEARS TO BE MUCH SLOWER THAN
% register.m
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 

%% basic variables
rr = o.ReferenceRound;
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';

%% loop through all tiles, finding spots in anchor channel on ref round
o.RawLocalYX = cell(nTiles,1);  % cell array, giving spots in local coordinates
o.RawIsolated = cell(nTiles,1);
SE = fspecial('disk', o.RegSmooth);

for t=NonemptyTiles(:)'
    [y,x] = ind2sub([nY nX], t);
    if mod(t,10)==0; fprintf('Loading tile %d anchor image\n', t); end
    
    AnchorIm  = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);
    if o.RegSmooth
        AnchorImSm = imfilter(AnchorIm ,SE);
    else
        AnchorImSm = AnchorIm;
    end
    [o.RawLocalYX{t}, o.RawIsolated{t}] = o.detect_spots(AnchorImSm);
end
%% get arrays ready


% WithinTileShift(t,:,r) is origin of tile t round r relative to origin of 
% tile t ref round
%WithinTileShift = nan(nTiles,2,o.nRounds);

% VerticalPairs: n x 2 array of tile IDs
% vShiftsYX: n x 2 array of YX shifts
% ccv: n x 1 array of correl coefs
% HorizontalPairs, hShiftsYX, cch: similar
VerticalPairs = zeros(0,2);
HorizontalPairs = zeros(0,2);
vShifts = zeros(0,2);
hShifts = zeros(0,2);
vScore = zeros(0,1);
hScore = zeros(0,1);
vChangedSearch = 0;
hChangedSearch = 0;

%% now do the alignments
for t=NonemptyTiles
    [y,x] = ind2sub([nY nX], t);    
    
    % can I align ref round to south neighbor?
    if y<nY && ~o.EmptyTiles(t+1)
        tic
        [shift, score, ChangedSearch] = o.get_initial_shift2(o.RawLocalYX{t}, o.RawLocalYX{t+1}, o.RegSearch.South,'Register');
        toc
        if all(isfinite(shift))
            VerticalPairs = [VerticalPairs; t, t+1];
            vShifts = [vShifts; shift];
            vScore = [vScore; score];
            vChangedSearch = vChangedSearch + ChangedSearch;
        end
        %ShowPos(o, y, x, y+1, x, rr, shift);
        fprintf('Tile %d (%d, %d), down: shift %d %d, score %f\n', t, y, x, shift, score);

        %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
        if size(vShifts,1) == 3 || (mod(vChangedSearch,2) == 0) && (vChangedSearch>0)
            o = o.GetNewSearchRange_Register(vShifts,'South');
        end
        
    end
    
    % can I align to east neighbor
    if x<nX && ~o.EmptyTiles(t+nY)
        tic
        [shift, score, ChangedSearch] = o.get_initial_shift2(o.RawLocalYX{t}, o.RawLocalYX{t+nY}, o.RegSearch.East,'Register');
        toc
        if all(isfinite(shift))
            HorizontalPairs = [HorizontalPairs; t, t+nY];
            hShifts = [hShifts; shift];
            hScore = [hScore; score];
            hChangedSearch = hChangedSearch + ChangedSearch;
        end        
        %ShowPos(o, y, x, y, x+1, rr, shift);
        fprintf('Tile %d (%d, %d), right: shift %d %d, score %f\n', t, y, x, shift, score);

        %Change search range after 3 tiles or if search has had to be widened twice (This is for speed).
        if size(hShifts,1) == 3 || (mod(hChangedSearch,2) == 0) && (hChangedSearch>0)
            o = o.GetNewSearchRange_Register(hShifts,'East');
        end
        
        
    end
               
end
tic
%Set any anomalous shifts to average of all other shifts
%Anomalous if awful score or either shift is an outlier
%Vertical
vOutlier = zeros(size(vShifts));
AnomalousScores = vScore<o.RegMinScore;
if max(AnomalousScores)>0
    warning('Looking at anomalous vShifts');
    AnomalousShift = max(isoutlier(vShifts(:,1)),isoutlier(vShifts(:,2)));
    AwfulScore = vScore < o.RegAbsoluteMinScore;
    for i=1:size(vScore,1)
        if min([AnomalousScores(i),AnomalousShift(i)+AwfulScore(i)])>0
            vOutlier(i,:) = vShifts(i,:);
            vShifts(i,:) = round(mean(vShifts(AnomalousScores==0,:)));    
            warning('vShift(%d) changed',i);
        end
    end
end

%Horizontal
hOutlier = zeros(size(hShifts));
AnomalousScores = hScore<o.RegMinScore;
if max(AnomalousScores)>0
    warning('Looking at anomalous hShifts');
    AnomalousShift = max(isoutlier(hShifts(:,1)),isoutlier(hShifts(:,2)));
    AwfulScore = hScore < o.RegAbsoluteMinScore;
    for i=1:size(hScore,1)
        if min([AnomalousScores(i),AnomalousShift(i)+AwfulScore(i)])>0
            hOutlier(i,:) = hShifts(i);
            hShifts(i,:) = round(mean(hShifts(AnomalousScores==0,:)));
            warning('hShift(%d) changed',i);
        end
    end
end
toc


%Save registration info for debugging
o.RegInfo.VerticalPairs = VerticalPairs;
o.RegInfo.vShifts = vShifts;
o.RegInfo.vScore = vScore;
o.RegInfo.vChangedSearch = vChangedSearch;
o.RegInfo.vOutlier = vOutlier;
o.RegInfo.HorizontalPairs = HorizontalPairs;
o.RegInfo.hShifts = hShifts;
o.RegInfo.hScore = hScore;
o.RegInfo.hChangedSearch = hChangedSearch;
o.RegInfo.hOutlier = hOutlier;
%save(fullfile(o.OutputDirectory, 'o2.mat'), 'o');

%% now we need to solve a set of linear equations for each shift,
% This will be of the form M*x = c, where x and c are both of length 
% nTiles=nY*nX. The t'th row is the equation for tile t. 
% c has columns for y and x coordinates

M = zeros(nTiles, nTiles);
c = zeros(nTiles, 2);
for i=1:size(VerticalPairs,1)
    if isnan(vShifts(i,1)); continue; end
    t1 = VerticalPairs(i,1);
    t2 = VerticalPairs(i,2);
    M(t1,t1) = M(t1,t1)+1;
    M(t1,t2) = M(t1,t2)-1;
    c(t1,:) = c(t1,:) - vShifts(i,:);
    M(t2,t2) = M(t2,t2)+1;
    M(t2,t1) = M(t2,t1)-1;
    c(t2,:) = c(t2,:) + vShifts(i,:);
end

for i=1:size(HorizontalPairs,1)
    if isnan(hShifts(i,1)); continue; end
    t1 = HorizontalPairs(i,1);
    t2 = HorizontalPairs(i,2);
    M(t1,t1) = M(t1,t1)+1;
    M(t1,t2) = M(t1,t2)-1;
    c(t1,:) = c(t1,:) - hShifts(i,:);
    M(t2,t2) = M(t2,t2)+1;
    M(t2,t1) = M(t2,t1)-1;
    c(t2,:) = c(t2,:) + hShifts(i,:);
end

% now we want to anchor one of the tiles to a fixed coordinate. We do this
% for a home tile in the middle, because it is going to be connected; and we set
% its coordinate to a large value, so any non-connected ones can be
% detected. (BTW this is why spectral clustering works!!)
Huge = 1e6;
TileDistFromCenter = abs(mod(0:nTiles-1, nY)-nY/2) + ...
    abs(floor((0:nTiles-1)/nY)-nX/2);
[~, HomeTile] = min(TileDistFromCenter(:)./~o.EmptyTiles(:));
%sub2ind([nY nX], ceil(nY/2), ceil(nX/2));
M(nTiles+1,HomeTile) = 1;
c(nTiles+1,:) = [Huge, Huge];

Tiny = 1e-4; % for regularization
TileOffset0 = (M+Tiny*eye(nTiles+1, nTiles))\c;

% find tiles that are connected to the home tile 
AlignedOK = (TileOffset0(:,1)>Huge/2);
TileOffset1 = nan(nTiles, 2);
TileOffset1(AlignedOK,:) = TileOffset0(AlignedOK,:)-Huge;

% RefPos(t,1:2) is origin of reference tile
RefPos = bsxfun(@minus,TileOffset1, nanmin(TileOffset1))+1;

% tile origin(t,1:2,r)
o.TileOrigin = zeros(nTiles,2,o.nRounds+o.nExtraRounds);
o.TileOrigin(:,:,rr) =  RefPos;

%%

%save(fullfile(o.OutputDirectory, 'o1.mat'), 'o');



%% now make background image
AnchorOrigin = round(o.TileOrigin(:,:,rr));
MaxTileLoc = max(AnchorOrigin);
BigDapiIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');
BigAnchorIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');

for t=NonemptyTiles
    MyOrigin = AnchorOrigin(t,:);
    if mod(t,10)==0; fprintf('Loading tile %d DAPI image\n', t); end
    if ~isfinite(MyOrigin(1)); continue; end
    LocalDapiIm = imread(o.TileFiles{o.ReferenceRound,t}, o.DapiChannel);
    BigDapiIm(floor(MyOrigin(1))+(1:o.TileSz), ...
        floor(MyOrigin(2))+(1:o.TileSz)) ...
        = imresize(LocalDapiIm, 1);
    LocalAnchorIm = imread(o.TileFiles{o.ReferenceRound,t}, o.AnchorChannel);
    BigAnchorIm(floor(MyOrigin(1))+(1:o.TileSz), ...
        floor(MyOrigin(2))+(1:o.TileSz)) ...
        = LocalAnchorIm;
end

o.BigDapiFile = fullfile(o.OutputDirectory, 'background_image.tif');

imwrite(BigDapiIm, o.BigDapiFile);
imwrite(BigAnchorIm, fullfile(o.OutputDirectory, 'anchor_image.tif'));

return
end
%%
function ShowPos(o, y, x, y1, x1, r, shift)
	if all(isfinite(shift))
		Color = 'b';
	else
		Color = 'r';
	end
    %figure(239856); 
    clf; hold on
    plot(o.TilePosYX(:,2), o.TilePosYX(:,1), 'k.');
    plot([x x1], [y y1], Color);
    plot(x, y, [Color 'o'], 'markersize', r*3);
    set(gca, 'ydir', 'reverse');
    title(sprintf('Round %d', r));
    drawnow;
end