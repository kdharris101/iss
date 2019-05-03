function o = call_spots(o)
% o = o.call_spots
% calls spots to codes for in-situ sequencing. Run this after find_spots
% 
% produces SpotGene{Spot}: name of gene for each spot
% SpotCode{Spot}: text representation of code for each spot 
% SpotScore(Spot): score saying how well the code fits (0...1)
% SpotIntensity(Spot): RMS intensity of the spot
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
  
nChans = o.nBP+2;

SpotColors = bsxfun(@rdivide, o.cSpotColors, prctile(o.cSpotColors, o.SpotNormPrctile));

% now we cluster the intensity vectors to estimate the Bleed Matrix
BleedMatrix = zeros(o.nBP,o.nBP,o.nRounds); % (Measured, Real, Round)
for r =1:o.nRounds
    m = squeeze(SpotColors(o.cSpotIsolated,:,r)); % data: nCodes by nBases
    
    [Cluster, v, s2] = ScaledKMeans(m, eye(4));
    for i=1:4
        BleedMatrix(:,i,r) = v(i,:) * sqrt(s2(i));
    end
end

if o.Graphics
    figure(98043765); clf
    for i=1:o.nRounds
        subplot(2,3,i); 
        imagesc(BleedMatrix(:,:,i)); 
        caxis([0 1]); 
        title(sprintf('Cycle %d', i)); 
        set(gca, 'xtick', 1:4);
        set(gca, 'XTickLabel', {'T', 'G', 'C', 'A'});
        set(gca, 'ytick', 1:4);
        set(gca, 'yTickLabel', {'T', 'G', 'C', 'A'});
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
    subplot(2,3,6);
    caxis([0 1]); 
    axis off
    colormap hot
    colorbar
    save(fullfile(o.OutputDirectory, 'BleedMatrix.mat'), 'BleedMatrix');
end

% now load in the code book and apply bleeds to it
codebook_raw = importdata(o.CodeFile);
CharCode = codebook_raw.textdata(2:end,5);
GeneName = codebook_raw.textdata(2:end,3);
nCodes = size(CharCode,1) - nnz(cellfun(@(v) strcmp(v(1:2),'SW'), CharCode)); % bit of a hack to get rid of Sst and Npy (assume always in the end)

% put them into object o but without the extras
o.CharCodes=CharCode(1:nCodes);
o.GeneNames=GeneName(1:nCodes);

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r = 1:o.nRounds
    if o.AnchorChannel == 2
        NumericalCode(:,r) = codebook_raw.data(1:nCodes,(r-1)*nChans + (o.AnchorChannel+1:nChans))*(1:o.nBP)';
    else
        NumericalCode(:,r) = codebook_raw.data(1:nCodes,(r-1)*nChans + (o.DapiChannel+1:nChans-1))*(1:o.nBP)';
    end
end

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:o.nRounds
        BledCodes(i,(1:o.nBP) + (r-1)*o.nBP) = BleedMatrix(:, NumericalCode(i,r), r);
        UnbledCodes(i,NumericalCode(i,r) + (r-1)*o.nBP) = 1;
    end
end

NormBledCodes = bsxfun(@rdivide, BledCodes, sqrt(sum(BledCodes.^2,2)));
FlatSpotColors = SpotColors(:,:);
o.SpotIntensity = sqrt(sum(FlatSpotColors.^2,2));
NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, o.SpotIntensity);
SpotScores = NormFlatSpotColors * NormBledCodes';

[o.SpotScore, BestCode] = max(SpotScores,[],2);
o.SpotCodeNo = uint16(BestCode);
o.SpotCombi = true(size(o.SpotCodeNo,1),1);

