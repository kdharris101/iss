function o = add_extra_genes(o)
% [GlobalYX, Gene] = iss_single_genes(o)
%
% processes additional rounds in which single genes were detected
% non-combinatorially. Outputs are spot positions in global coordinates 
% (see iss_find_spots) and Gene assignments (see iss_call_spots)
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 
%% basic variables
rr = o.ReferenceRound;
Tiles = find(~o.EmptyTiles)';

[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
nExtraCodes = size(o.ExtraCodes,1);
nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));

%% first put the gene names into the structure
for i=1:size(o.ExtraCodes,1)
    o.CharCodes{nCombiCodes+i} = 'EXTRA';
    o.GeneNames{nCombiCodes+i} = o.ExtraCodes{i,1};
end

%% main loop
for i=1:nExtraCodes
    for t=Tiles(:)'
        if mod(t,10)==0; fprintf('extra %s spots for tile %d\n', o.ExtraCodes{i,1}, t); end
        r = o.ExtraCodes{i,2};
        c = o.ExtraCodes{i,3};
        ExtraIm = imread(o.TileFiles{r,t}, c); 
        ExtraRawLocalYX{t,i} = o.detect_spots(ExtraIm);
        ExtraRawGlobalYX{t,i} = bsxfun(@plus, ExtraRawLocalYX{t,i}, o.TileOrigin(t,:,r));
        nSpots = size(ExtraRawGlobalYX{t,i},1);
        ExtraRawCodeNo{t,i} = repmat(nCombiCodes+i, nSpots,1);
        ExtraRawDetectedTile{t,i} = repmat(t,nSpots,1);
        
        if o.SmoothSize
            ExtraImSm = imfilter(double(ExtraIm), fspecial('disk', o.SmoothSize));
        else
            ExtraImSm = BaseIm;
        end
        ExtraRawIntensity{t,i} = IndexArrayNan(ExtraImSm, ExtraRawLocalYX{t,i}');
    end
end

%% concatenate big arrays
ExtraAllGlobalYX = vertcat(ExtraRawGlobalYX{:});
ExtraAllDetectedTile = vertcat(ExtraRawDetectedTile{:});
ExtraAllCodeNo = vertcat(ExtraRawCodeNo{:});
ExtraAllIntensity = vertcat(ExtraRawIntensity{:});

%% normalize intensities
% ExtraAllIntensityNorm = zeros(size(ExtraAllIntensity));
% for i=1:nExtraCodes
%     MySpots = (ExtraAllCodeNo==nCombiCodes+i);
%     DivideBy = prctile(ExtraAllIntensity(MySpots), o.SpotNormPrctile);
%     ExtraAllIntensityNorm(MySpots) = ExtraAllIntensity(MySpots)/DivideBy;
% end
  

%% eliminate duplicates (can go straight from all to good here, no nd)
% note that by this stage everything is in global coordinates relative to
% reference frame, which is why we don't need to add to o.RefPos
ExtraAllTile = which_tile(ExtraAllGlobalYX, o.TileOrigin(:,:,rr), o.TileSz);
ExtraGood = (ExtraAllTile==ExtraAllDetectedTile);



%% now writing the cells, going after the last combi cell
nCombiCells = sum(o.SpotCombi);
NewCellIndices = (nCombiCells+1):(nCombiCells + sum(ExtraGood));

o.SpotGlobalYX(NewCellIndices,:) =ExtraAllGlobalYX(ExtraGood,:);
o.SpotCodeNo(NewCellIndices) = ExtraAllCodeNo(ExtraGood);
o.SpotCombi(NewCellIndices) = false;
o.SpotScore(NewCellIndices) = 1;
o.SpotIntensity(NewCellIndices) = ExtraAllIntensity(ExtraGood);

