function o=register(o)
% o=iss_register(o)
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
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
%% get arrays ready
rr = o.ReferenceRound;

% we index tiles by their xy coordinate (sometimes as a linear index). Not
% all of these tiles are actually there. NonemptyTiles lists the ones that
% are.
NonemptyTiles = find(~o.EmptyTiles)';

[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;

% stores pre-computed FFT and energy (to save time)
RefFFTStore = cell{nTiles,1};


%% first we align rounds on each tile individually
% WithinTileShift(r,:,t) is origin of tile t round r relative to origin of 
% tile t ref round
WithinTileShift = nan(nTiles,o.nRounds);

% linear shift
for t=NonemptyTiles
    [y,x] = ind2sub([nY nX], t);
    
    if mod(t,10)==0; fprintf('Loading tile %d anchor image\n', t); end
    RefImage = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);

    for r=1:o.nRounds+o.nExtraRounds
        
        if r==rr % no offset for reference round
            WithinTileShift(r, :, t) = [0 0];
            continue;
        end
 
        MyTile = imread(o.TileFiles{r,y,x},o.AnchorChannel);

        % align to same tile in reference round, cacheing fft if not
        % already there
        if isempty(RefFFTStore{t})
            [shift, cc, f] = ImRegFft2(RefImage, MyTile, o.RegCorrThresh, o.RegMinSize);
            RefFFTStore{t} = f;
        else
            [shift, cc] = ImRegFft2(RefFFTStore{t}, MyTile, o.RegCorrThresh, o.RegMinSize);
        end
        ShowPos(o, y, x, y, x, r, shift);
        fprintf('\nround %d, tile %d at (%d, %d): shift %d %d, to ref round, cc %f\n', r, t, y, x, shift, cc);
        WithinTileShift(r,:,t) = shift;
     
    
    end
    save o2 o
end


%% now we stitch neighboring tiles on the reference round

% first make a list of vertical and horizontal tile pairs
dy =1; dx = nY;

PotentialVerticalPairs = [1:nT; (1:nT)+dy]';
BothThere = all(ismember(PotentialVerticalPairs,NonemptyTiles),2);
NotBot = ones(nT,1); NotBot(nY:nY:end)=0;
VerticalPairs = PotentialVerticalPairs(BothThere&NotBot,:);

PotentialHorizontalPairs = [1:nT; (1:nT)+dx]';
BothThere = all(ismember(PotentialHorizontalPairs,NonemptyTiles),2);
NotRt = ones(nT,1); NotRt(nY*(nX-1)+1:end)=0;
HorizontalPairs = PotentialHorizontalPairs(BothThere&NotRt,:);

nVerticalPairs = size(VerticalPairs,1);
nHorizontalPairs = size(HorizontalPairs,1);

% to store results of pairwise image registrations
% stores global coordinate of lower or right tile relative to upper or left
vShifts = nan(nVerticalPairs,2);
hShifts = nan(nHorizontalPairs,2);
ccv = zeros(nVerticalPairs,1);
cch = zeros(nHorizontalPairs,1);

for i=1:size(VerticalPairs,1)
    [y,x] = ind2sub([nY nX], VerticalPairs(i,1));
    [vShifts(i,:), ccv(i)] = ImRegFft2(RefFFTStore{VerticalPairs(i,1)}...
        ,RefFFTStore{VerticalPairs(i,1)}, o.RegCorrThresh, o.RegMinSize);
    ShowPos(o, y, x, y+1, x, rr, vShifts(i,:));
    fprintf('%d, %d, down: shift %d %d, cc %f\n', y, x, vShifts(i,:), ccv(i));
end

for i=1:size(HorizontalPairs,1)
    [y,x] = ind2sub([nY nX], HorizontalPairs(i,1));
    [hShifts(i,:), cch(i)] = ImRegFft2(RefFFTStore{HorizontalPairs(i,1)}...
        ,RefFFTStore{HorizontalPairs(i,1)}, o.RegCorrThresh, o.RegMinSize);
    ShowPos(o, y, x, y+1, x, rr, vShifts(i,:));
    fprintf('%d, %d, down: shift %d %d, cc %f\n', y, x, vShifts(i,:), ccv(i));
end

% THAT'S AS FAR AS I GOT THIS MORNING ...

% now we need to solve a set of linear equations for each shift, that makes
% each tile position the mean of what is suggested by all pairs it appears
% in. This will be of the form M*x = c, where x and c are both of length 
% nTiles=nY*nX. The t'th row is the equation for tile t. 
% c has columns for y and x coordinates

M = zeros(nTiles, nTiles);
c = zeros(nTiles, 2);
for i=find(VerticalPairs)'
    if isnan(vShifts(i,1)); continue; end
    [y1,x1] = ind2sub(size(VerticalPairs), i);
    y2 = y1+1; x2 = x1;
    t1 = sub2ind([nY nX], y1, x1);
    t2 = sub2ind([nY nX], y2, x2);
    M(t1,t1) = M(t1,t1)+1;
    M(t1,t2) = M(t1,t2)-1;
    c(t1,:) = c(t1,:) - vShifts(i,:);
    M(t2,t2) = M(t2,t2)+1;
    M(t2,t1) = M(t2,t1)-1;
    c(t2,:) = c(t2,:) + vShifts(i,:);
end

for i=find(HorizontalPairs)'
    if isnan(hShifts(i,1)); continue; end
    [y1,x1] = ind2sub(size(HorizontalPairs), i);
    y2 = y1; x2 = x1+1;
    t1 = sub2ind([nY nX], y1, x1);
    t2 = sub2ind([nY nX], y2, x2);
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
o.RefPos = bsxfun(@minus,TileOffset1, nanmin(TileOffset1))+1;

save o1 o




%% now make background image

MaxTileLoc = max(o.RefPos);
BigDapiIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');
BigAnchorIm = zeros(ceil((MaxTileLoc + o.TileSz)), 'uint16');

for t=NonemptyTiles
    if mod(t,10)==0; fprintf('Loading tile %d anchor image\n', t); end
    if ~isfinite(o.RefPos(t,1)); continue; end
    LocalDapiIm = imread(o.TileFiles{o.ReferenceRound,t}, o.DapiChannel);
    BigDapiIm(floor(o.RefPos(t,1))+(1:o.TileSz), ...
        floor(o.RefPos(t,2))+[1:o.TileSz]) ...
        = imresize(LocalDapiIm, 1);
    LocalAnchorIm = imread(o.TileFiles{o.ReferenceRound,t}, o.AnchorChannel);
    BigAnchorIm(floor(o.RefPos(t,1))+(1:o.TileSz), ...
        floor(o.RefPos(t,2))+(1:o.TileSz)) ...
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
    plot(x, y, [Color 'o']);
    set(gca, 'ydir', 'reverse');
    title(sprintf('Round %d', r));
    drawnow;
end