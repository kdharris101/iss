function [o, VerticalPairs, vShifts, HorizontalPairs, hShifts] = get_Fft_shifts(o)
%% [o, VerticalPairs, vShifts, HorizontalPairs, hShifts] = get_Fft_shifts(o)
%
% function for finding shifts between overlapping tiles in registration
% step using Fft method.

%% load and store ref images
% we index tiles by their yx coordinate (sometimes as a linear index). Not
% all of these tiles are actually there. NonemptyTiles lists the ones that
% are.

rr = o.ReferenceRound;
[nY, nX] = size(o.EmptyTiles);
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end
RefImages = zeros(o.TileSz, o.TileSz, nY, nX, 'uint16');

for t=NonemptyTiles(:)'
    [y,x] = ind2sub([nY nX], t);
    if mod(t,10)==0; fprintf('Loading tile %d anchor image\n', t); end
    Im = imread(o.TileFiles{rr,y,x}, o.ReferenceChannel);
    if o.RegSmooth
        RefImages(:,:,t) = imfilter(Im, fspecial('disk', o.RegSmooth));
    else
        RefImages(:,:,t) = Im;
    end
end

%% get arrays ready


% VerticalPairs: n x 2 array of tile IDs
% vShiftsYX: n x 2 array of YX shifts
% vScore: n x 1 array of correl coefs
% HorizontalPairs, hShiftsYX, cch: similar

VerticalPairs = zeros(0,2);
HorizontalPairs = zeros(0,2);
vShifts = zeros(0,2);
hShifts = zeros(0,2);
vScore = zeros(0,1);
hScore = zeros(0,1);
Graphics = o.Graphics;
BadShifts = 0;
AbsoluteMaxShift = o.RegInfo.AbsoluteMaxShift;
AbsoluteMinShift = o.RegInfo.AbsoluteMinShift;

%% now do the alignments
for t=NonemptyTiles
    if t == NonemptyTiles(1)
        o.Graphics = 2;     
    end
    [y,x] = ind2sub([nY nX], t);    
    
    % can I align ref round to south neighbor?
    if y<nY && ~o.EmptyTiles(t+1)
        tic
        %Works best if get rid of negative values i.e. unit16(-10) = 0 as done here
        [shift, cc] = o.ImRegFft2_Register(RefImages(:,:,t)-o.TilePixelValueShift,...
            RefImages(:,:,t+1)-o.TilePixelValueShift, o.RegCorrThresh, o.RegMinSize,'South');
        if all(isfinite(shift))
            VerticalPairs = [VerticalPairs; t, t+1];
            vShifts = [vShifts; shift];
            vScore = [vScore; cc];
        end
        %ShowPos(o, y, x, y+1, x, rr, shift);
        fprintf('Tile %d (%d, %d), down: shift %d %d, cc %f\n', t, y, x, shift, cc);
        toc
    end
    
    % can I align to east neighbor
    if x<nX && ~o.EmptyTiles(t+nY)
        tic
        [shift, cc] = o.ImRegFft2_Register(RefImages(:,:,t)-o.TilePixelValueShift,...
            RefImages(:,:,t+nY)-o.TilePixelValueShift, o.RegCorrThresh, o.RegMinSize,'East');
        if all(isfinite(shift))
            HorizontalPairs = [HorizontalPairs; t, t+nY];
            hShifts = [hShifts; shift];
            hScore = [hScore; cc];
        end        
        %ShowPos(o, y, x, y, x+1, rr, shift);
        fprintf('Tile %d (%d, %d), right: shift %d %d, cc %f\n', t, y, x, shift, cc);
        toc
    end
            
    if t == NonemptyTiles(1)
        o.Graphics = Graphics;
    end
    
end

o.RegInfo.Fft.VerticalPairs = VerticalPairs;
o.RegInfo.Fft.HorizontalPairs = HorizontalPairs;
o.RegInfo.Fft.vShifts = vShifts;
o.RegInfo.Fft.hShifts = hShifts;
o.RegInfo.Fft.vScore = vScore;
o.RegInfo.Fft.hScore = hScore;