function o = extract_and_filter_parallel(o)
%% This part is the parallel stuff
q = parallel.pool.DataQueue;
nTotalRounds = o.nRounds+o.nExtraRounds;
parpool('local',nTotalRounds);
p = gcp();

for idx = 1:nTotalRounds
  f(idx) = parfeval(p,@get_tile_pos,2,o,idx,q); 
end

TilePosYX = cell(1,nTotalRounds);
fName = cell(1,nTotalRounds);
for idx = 1:nTotalRounds
  % fetchNext blocks until next results are available.
  [completedIdx,value1,value2] = fetchNext(f);
  listener = q.afterEach(@disp);            %NOT QUITE SURE HOW THIS DISPLAYING IS DONE
  TilePosYX{completedIdx} = value1;
  fName{completedIdx} = value2;
  fprintf('Finished Round %d.\n', completedIdx);
end

%Clean up all data outputted from parallel loops
o.TilePosYX = TilePosYX{1};
nTiles = size(o.TilePosYX,1);
o.TileFiles = cell(o.nRounds+o.nExtraRounds,max(o.TilePosYX(:,1)),max(o.TilePosYX(:,2)));
for r=1:nTotalRounds
    for t=1:nTiles
        o.TileFiles{r,o.TilePosYX(t,1), o.TilePosYX(t,2)} = fName{r}{t};
    end
end

o.EmptyTiles = cellfun(@isempty, squeeze(o.TileFiles(o.ReferenceRound,:,:)));

% a = get_tile_pos(o,8);
end
%% This part is the same as extract_and_filter.m
function [TilePosYX,fName] = get_tile_pos(o,r,q)

imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);

% construct a Bio-Formats reader with the Memoizer wrapper
bfreader = loci.formats.Memoizer(bfGetReader(), 0);
% initiate reader
bfreader.setId(imfile);

% get some basic image metadata
[nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
    get_ome_tilepos(bfreader);
if isempty(xypos) || size(xypos, 1)==1
    if r == 1
        warning('first round xypos empty - using values from initial manual input')
        assert(~isempty(o.TileInitialPosXY), 'xypos unavailable')
        xypos = o.TileInitialPosXY;
        xyposOld = xypos;
    else
        warning('xypos empty - using values from previous round')
        xypos = xyposOld;
    end
    nSerieswPos = size(xypos,1);
else
    xyposOld = xypos;
end

scene = nSeries/nSerieswPos;

bfreader.close();


% find x and y grid spacing as median of distances that are about
% right
dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
dy = xypos(:,1)-xypos(:,1)'; % all pairs of y distances
yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));


% find coordinates for each tile
o.TileInitialPosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
TilePosYX = o.TileInitialPosYX;
%Below is a safeguard incase wrong positions found - can do
%this as we know what the answer should be.
MaxY = max(TilePosYX(:,1));
MaxX = max(TilePosYX(:,2));
if MaxY*MaxX ~= nSeries
    warning('Number of tiles (%d) is not equal to maximum Y position (%d) multiplied by maximum X position (%d)'...
        , nSeries, MaxY, MaxX)
else
    TilePosY = flip(repelem(1:MaxY,MaxX));
    TilePosYX(:,1) = TilePosY;
    TilePosX = repmat([flip(1:MaxX),1:MaxX],1,ceil(MaxY/2));
    TilePosYX(1:nSeries,2) = TilePosX(1:nSeries);
end
%send(q,fprintf('Round %d tile %d colour channel %d finished.\n', r, 1, 1));


% set up filename grid for this round
fName = cell(nSerieswPos,1);

if strcmpi(o.ExtractR,'auto')
    o.ExtractR = round(1/pixelsize);
end
if strcmpi(o.DapiR,'auto')
    o.DapiR = round(8/pixelsize);
end
SE = strel('disk', o.ExtractR);
DapiSE = strel('disk', o.DapiR);

for t=1:nSerieswPos 
    
    fName{t} = fullfile(o.TileDirectory, ...
        [o.FileBase{r}, '_t', num2str(t), '.tif']);
    
    if exist(fName{t}, 'file')       
        q.send('Round '+string(r)+' tile '+string(t)+' already done');        
        continue;
    end
    
    bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
    % use the memo file cached before
    bfreader.setId(imfile);
    bfreader.setSeries(scene*t-1);
    
    for c=1:nChannels
        %Only need anchor and dapi tiles in reference round
        if (r == o.ReferenceRound && c ~= o.AnchorChannel) && (r == o.ReferenceRound && c ~= o.DapiChannel)
            IFS = zeros(o.TileSz,o.TileSz);
        else
            I = cell(nZstacks,1);
            for z = 1:nZstacks
                iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                I{z} = bfGetPlane(bfreader, iPlane);
            end
            
            % focus stacking
            IFS = o.fstack_modified(I);
            
            % tophat
            if c == o.DapiChannel && r == o.ReferenceRound
                IFS = imtophat(IFS, DapiSE);
            else
                IFS = imtophat(IFS, SE);
            end
        end
        % write stack image
        imwrite(uint16(IFS),...
            fullfile(o.TileDirectory,...
            [o.FileBase{r}, '_t', num2str(t), '.tif']),...
            'tiff', 'writemode', 'append');
    end
    q.send('Round '+string(r)+' tile '+string(t)+' finished'); 
    bfreader.close();
end

end
