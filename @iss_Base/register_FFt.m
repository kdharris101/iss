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
 

%% load and store ref images
% we index tiles by their yx coordinate (sometimes as a linear index). Not
% all of these tiles are actually there. NonemptyTiles lists the ones that
% are.

rr = o.ReferenceRound;
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
RefImages = zeros(o.TileSz, o.TileSz, nY, nX, 'uint16');

for t=NonemptyTiles(:)'
    [y,x] = ind2sub([nY nX], t);
    if mod(t,10)==0; fprintf('Loading tile %d anchor image\n', t); end
    Im = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);
    if o.RegSmooth
        RefImages(:,:,t) = imfilter(Im, fspecial('disk', o.RegSmooth));
    else
        RefImages(:,:,t) = Im;
    end
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
ccv = zeros(0,1);
cch = zeros(0,1);

%% now do the alignments
for t=NonemptyTiles
    [y,x] = ind2sub([nY nX], t);    
    
    % can I align ref round to south neighbor?
    if y<nY && ~o.EmptyTiles(t+1)
        tic
        [shift, cc] = o.ImRegFft2_Register(RefImages(:,:,t), RefImages(:,:,t+1), o.RegCorrThresh, o.RegMinSize,'South');
        if all(isfinite(shift))
            VerticalPairs = [VerticalPairs; t, t+1];
            vShifts = [vShifts; shift];
            ccv = [ccv; cc];
        end
        %ShowPos(o, y, x, y+1, x, rr, shift);
        fprintf('Tile %d (%d, %d), down: shift %d %d, cc %f\n', t, y, x, shift, cc);
        toc
    end
    
    % can I align to east neighbor
    if x<nX && ~o.EmptyTiles(t+nY)
        tic
        [shift, cc] = o.ImRegFft2_Register(RefImages(:,:,t), RefImages(:,:,t+nY), o.RegCorrThresh, o.RegMinSize,'East');
        if all(isfinite(shift))
            HorizontalPairs = [HorizontalPairs; t, t+nY];
            hShifts = [hShifts; shift];
            cch = [cch; cc];
        end        
        %ShowPos(o, y, x, y, x+1, rr, shift);
        fprintf('Tile %d (%d, %d), right: shift %d %d, cc %f\n', t, y, x, shift, cc);
        toc
    end
            
    
    
end

o.RegInfo.VerticalPairs = VerticalPairs;
o.RegInfo.HorizontalPairs = HorizontalPairs;
o.RegInfo.vShifts = vShifts;
o.RegInfo.hShifts = hShifts;
o.RegInfo.ccv = ccv;
o.RegInfo.cch = cch;

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

o.BigDapiFile = fullfile(o.OutputDirectory, 'background_image.tif');
imwrite(BigDapiIm, o.BigDapiFile);
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