function [SpotColors,PointCorrectedLocalYX] = get_spot_colors(o,LocalYX,LocalTile,...
    RoundTile,AllBaseLocalYX,nMatches,Error)
%% GET_SPOT_COLORS Summary of this function goes here
% loop through all tiles, finding spot colors
% o: iss object
% LocalYX(s): YX coordinate of spot s on tile given by LocalTile(s).
% LocalTile(s): Tile that spot s was found on.
% RoundTile(s,r): stores appropriate tile for spot s on round r
% AllBaseLocalYX: cell(nTiles,o.nBP,nRounds). AllBaseLocalYX{t,b,r} gives
% YX local coordinates of all spots found on tile t, channel b, round r.
% nMatches(t,b,r): gives PCR matches to tile t, channel b, round r.
% Error(t,b,r): gives PCR error to tile t, channel b, round r. 
% SpotColors(s,b,r): gives intensity for spot s in channel b, round r.
% PointCorrectedLocalYX(s,:,r,b): gives local coordinate of spot s but
% transformed to round r, channel b according to o.D(:,:,t,r,b).

Verbose = nargin>=6;        %print if have more than 6 arguments. 

NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end
[nY, nX] = size(o.EmptyTiles);
AllBaseSpotNo = cell2mat(cellfun(@size,AllBaseLocalYX,'uni',false));
AllBaseSpotNo = AllBaseSpotNo(:,1:2:o.nBP*2,:);

nSpots = size(LocalYX,1);
LocalYX = [LocalYX-o.TileCentre,ones(nSpots,1)];
SpotColors = nan(nSpots, o.nBP, max(o.UseRounds));
PointCorrectedLocalYX = nan(nSpots, 2, o.nRounds, o.nBP);

for t=NonemptyTiles
    [y, x] = ind2sub([nY nX], t);
   
    for r=o.UseRounds         
        % find spots whose home tile on round r is t
        MySpots = (RoundTile(:,r)==t);
        if ~any(MySpots); continue; end
        
        % open file for this tile/round
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = LocalTile(RoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        if Verbose
            fprintf('\nRef round home tiles for spots in t%d at (%2d, %2d), r%d: ', t, y, x, r);
            for i=MyRefTiles(:)'
                fprintf('t%d, %d spots; ', i, sum(RefRoundHomeTiles==i));
            end
            fprintf('\n');
        end
                
        % now read in images for each base
        for b=o.UseChannels                      
            TifObj.setDirectory(o.FirstBaseChannel + b - 1);
            BaseIm = int32(TifObj.read())-o.TilePixelValueShift;            
            if o.SmoothSize
                BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
            else
                BaseImSm = BaseIm;
            end
            
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (RoundTile(:,r)==t & LocalTile==t2);
                MyLocalYX = LocalYX(MyBaseSpots,:);
                
                if t == t2
                    if Verbose
                        fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                            t, t2, r, b,  nMatches(t,b,r), o.RawLocalNo(t2), Error(t,b,r));
                        if nMatches(t,b,r)<o.MinPCMatchFract*AllBaseSpotNo(t,b,r) || isempty(nMatches(t,b,r))
                            warning('Tile %d, channel %d, round %d has %d point cloud matches, which is below the threshold of %d.',...
                                t,b,r,nMatches(t,b,r),o.MinPCMatchFract*AllBaseSpotNo(t,b,r));
                        end
                    end
                    MyPointCorrectedYX = MyLocalYX*o.D(:,:,t,r,b)+o.TileCentre;
                    MyPointCorrectedYX = round(MyPointCorrectedYX);
                    PointCorrectedLocalYX(MyBaseSpots,:,r,b) = MyPointCorrectedYX;
                    SpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                else
                    [MyPointCorrectedYX, Error_diff_tile, nMatches_diff_tile] = o.different_tile_transform(AllBaseLocalYX,o.RawLocalYX, ...
                        MyLocalYX,t,t2,r,b);
                    if Verbose
                        fprintf('Point cloud: ref round tile %d -> tile %d round %d base %d, %d/%d matches, error %f\n', ...
                            t, t2, r, b,  nMatches_diff_tile, o.RawLocalNo(t2), Error_diff_tile);
                    end
                    if isempty(nMatches_diff_tile) || nMatches_diff_tile<o.MinPCMatchFract*AllBaseSpotNo(t,b,r)
                        continue;
                    end
                    PointCorrectedLocalYX(MyBaseSpots,:,r,b) = MyPointCorrectedYX;
                    SpotColors(MyBaseSpots,b,r) = IndexArrayNan(BaseImSm, MyPointCorrectedYX');
                end
               
            end    
        end
        TifObj.close();       
    end
end
if Verbose
    fprintf('\n');
end

end

