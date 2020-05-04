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
 
%% Logging

if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end

%% basic variables
rr = o.ReferenceRound;
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';

%Specify which rounds/colour channels to use (default is all)
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end

if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

if isempty(o.ReferenceSpotChannels)
    if rr == o.AnchorRound
        o.ReferenceSpotChannels = o.AnchorChannel;
    else
        o.ReferenceSpotChannels = o.UseChannels;
    end
end

%% loop through all tiles, finding spots in ref round
%Need spots in all channels as in theory, each spot should only appear in
%one channel
o.RawLocalYX = cell(nTiles,o.nBP);  % cell array, giving spots in local coordinates
o.RawIsolated = cell(nTiles,o.nBP);
o.RawLocalNo = zeros(nTiles,1);

for t=NonemptyTiles(:)'
    if mod(t,10)==0; fprintf('Loading tile %d reference image\n', t); end
    FileName = o.TileFiles{rr,t};
    TifObj = Tiff(FileName);
    for b=o.ReferenceSpotChannels        
        TifObj.setDirectory(o.FirstBaseChannel + b - 1);
        ReferenceIm = int32(TifObj.read())-o.TilePixelValueShift;            
        if o.SmoothSize
            SE = fspecial('disk', o.SmoothSize);
            ReferenceImSm = imfilter(ReferenceIm ,SE);
        else
            ReferenceImSm = ReferenceIm;
        end
        [o.RawLocalYX{t,b}, o.RawIsolated{t,b}] = o.detect_spots(ReferenceImSm,t,b,rr);
    end
    o.RawLocalNo(t) = length(vertcat(o.RawIsolated{t,:}));
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
AbsoluteMaxShift = max(o.RegSearch.South.Y);
AbsoluteMinShift = min(o.RegSearch.South.Y);
BadShifts = 0;
InSignificantTilesLeft = o.RegInsignificantFract*length(NonemptyTiles);
o.RegMaxBadShifts = min(o.RegMaxBadShifts,InSignificantTilesLeft);
o.RegInfo.SingleFft.Pairs = [];
o.RegInfo.SingleFft.OldShifts = [];
o.RegInfo.SingleFft.OldScores = [];

%% now do the alignments
for t=NonemptyTiles
    if isprop(o,'RegMethod') && strcmpi(o.RegMethod, 'Fft')
        break;
    end
    [y,x] = ind2sub([nY nX], t);    
    
    % can I align ref round to south neighbor?
    if y<nY && ~o.EmptyTiles(t+1)
        tic
        [shift, score, ChangedSearch] = o.get_initial_shift2(o.RawLocalYX{t,o.ReferenceChannel},...
            o.RawLocalYX{t+1,o.ReferenceChannel}, o.RegSearch.South,'Register');
        toc
        
        if shift(1) >= AbsoluteMaxShift || shift(1) <= AbsoluteMinShift ||...
                min(size(o.RawLocalYX{t,o.ReferenceChannel},1),...
                size(o.RawLocalYX{t+1,o.ReferenceChannel},1))<o.OutlierMinScore
            BadShifts = BadShifts+1;
            if BadShifts>=o.RegMaxBadShifts &&...
                    length(NonemptyTiles)-(find(NonemptyTiles==t)-1)>InSignificantTilesLeft
                %Only break if significant amount of tiles left
                break;
            end
            warning('tile %d to tile %d shift = [%d, %d], which is faulty, trying with Fft method',...
                t, t+1, shift(1), shift(2));
            o.RegInfo.SingleFft.Pairs = [o.RegInfo.SingleFft.Pairs; t, t+1];
            o.RegInfo.SingleFft.OldShifts = [o.RegInfo.SingleFft.OldShifts; shift(1), shift(2)];
            o.RegInfo.SingleFft.OldScores = [o.RegInfo.SingleFft.OldScores; score];
            [shift, score] = o.get_Fft_shift_single(t,o.ReferenceRound,o.ReferenceChannel,...
                t+1,o.ReferenceRound,o.ReferenceChannel,'Register','South');
        end
            
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
        [shift, score, ChangedSearch] = o.get_initial_shift2(o.RawLocalYX{t,o.ReferenceChannel},...
            o.RawLocalYX{t+nY,o.ReferenceChannel}, o.RegSearch.East,'Register');
        toc
        
        if shift(2) >= AbsoluteMaxShift || shift(2) <= AbsoluteMinShift ||...
                min(size(o.RawLocalYX{t,o.ReferenceChannel},1),...
                size(o.RawLocalYX{t+nY,o.ReferenceChannel},1))<o.OutlierMinScore
            BadShifts = BadShifts+1;
            if BadShifts>=o.RegMaxBadShifts &&...
                    length(NonemptyTiles)-(find(NonemptyTiles==t)-1)>InSignificantTilesLeft
                %Only break if significant amount of tiles left
                break;
            end
            warning('tile %d to tile %d shift = [%d, %d], which is faulty, trying with Fft method',...
                t, t+nY, shift(1), shift(2));
            o.RegInfo.SingleFft.Pairs = [o.RegInfo.SingleFft.Pairs; t, t+nY];
            o.RegInfo.SingleFft.OldShifts = [o.RegInfo.SingleFft.OldShifts; shift(1), shift(2)];
            o.RegInfo.SingleFft.OldScores = [o.RegInfo.SingleFft.OldScores; score];
            [shift, score] = o.get_Fft_shift_single(t,o.ReferenceRound,o.ReferenceChannel,...
                t+nY,o.ReferenceRound,o.ReferenceChannel,'Register','East');
        end
 
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

if t == NonemptyTiles(end)
    %Set any anomalous shifts to average of all other shifts
    %Anomalous if awful score or either shift is an outlier
    [vShifts, vOutlier] = o.AmendShifts(vShifts,vScore,'Register');
    [hShifts, hOutlier] = o.AmendShifts(hShifts,hScore,'Register');
    o.RegInfo.vOutlier = vOutlier;
    o.RegInfo.hOutlier = hOutlier;
end    
    
%Save registration info for debugging
o.RegInfo.Method = 'PointBased';
o.RegInfo.VerticalPairs = VerticalPairs;
o.RegInfo.vShifts = vShifts;
o.RegInfo.vScore = vScore;
o.RegInfo.vChangedSearch = vChangedSearch;
o.RegInfo.HorizontalPairs = HorizontalPairs;
o.RegInfo.hShifts = hShifts;
o.RegInfo.hScore = hScore;
o.RegInfo.hChangedSearch = hChangedSearch;
o.RegInfo.AbsoluteMaxShift = AbsoluteMaxShift;
o.RegInfo.AbsoluteMinShift = AbsoluteMinShift;
o.RegInfo.BadShifts = BadShifts;
%save(fullfile(o.OutputDirectory, 'o2.mat'), 'o');

if t ~= NonemptyTiles(end)
    warning(['The threshold number of %d tiles with shifts outside the range'...
        '%d to %d has been reached. Starting again with the Fft method.'],...
        BadShifts,AbsoluteMinShift,AbsoluteMaxShift)
    %If had to many bad shifts, try again with Fft method
    o.RegInfo.Method = 'Fft';
    [o, VerticalPairs, vShifts, HorizontalPairs, hShifts] = get_Fft_shifts(o);
end

if strcmpi(o.RegInfo.Method, 'Fft')
    %If point based method failed here, it will probably fail for
    %find_spots as well so use Fft method there.
    o.FindSpotsMethod = o.RegInfo.Method;   
end

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

%Plot stitched anchor image
if o.Graphics
    SquareX1 = [0, 0, o.TileSz];
    SquareY1 = [o.TileSz, 0, 0];
    SquareX2 = [o.TileSz, o.TileSz, 0];
    SquareY2 = [0, o.TileSz, o.TileSz];
    figure(53277);imagesc(BigAnchorIm-o.TilePixelValueShift);
    hold on
    for t=NonemptyTiles
        MyOrigin = o.TileOrigin(t,:,o.ReferenceRound);
        plot(SquareX1 + MyOrigin(2), SquareY1 + MyOrigin(1),...
            '--', 'Color', 'w');
        plot(SquareX2 + MyOrigin(2), SquareY2 + MyOrigin(1),...
            ':', 'Color', 'w');
        text(MyOrigin(2), MyOrigin(1),...
            sprintf('T%d r%d c%d', t, o.ReferenceRound, o.ReferenceChannel), 'color', 'w');
    end
    title('Stitched Anchor Image');
    hold off
end


if o.ReferenceRound == o.AnchorRound
    %Plot DapiImage showing tiling
    if o.Graphics
        figure(53278);imagesc(BigDapiIm);
        hold on
        for t=NonemptyTiles
            MyOrigin = o.TileOrigin(t,:,o.ReferenceRound);
            plot(SquareX1 + MyOrigin(2), SquareY1 + MyOrigin(1),...
                '--', 'Color', 'w');
            plot(SquareX2 + MyOrigin(2), SquareY2 + MyOrigin(1),...
                ':', 'Color', 'w');
            text(MyOrigin(2), MyOrigin(1),...
                sprintf('T%d r%d c%d', t, o.ReferenceRound, o.DapiChannel), 'color', 'w');
        end
        title('Stitched Dapi Image');
        hold off
    end
    o.BigDapiFile = fullfile(o.OutputDirectory, 'background_image.tif');
    imwrite(BigDapiIm, o.BigDapiFile);
end
o.BigAnchorFile = fullfile(o.OutputDirectory, 'anchor_image.tif');
imwrite(BigAnchorIm, o.BigAnchorFile);

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