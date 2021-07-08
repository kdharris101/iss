function o = get_omp_bled_codes(o)
%%  o = o.get_omp_bled_codes;
%   This gets bleed matrix for z-scored spot colors.
%   Then finds background codes, from pixel colors that can't be explained
%   by any genes. 
%   ompBledCodes is then the normalised z-scored bled codes for the genes
%   plus the best few background eigenvectors. 
%   o.Max_nBackground: there will be at most this number of background codes.
%   o.BackgroundMaxGeneDotProductThresh: will reject background codes that
%   have a dot product with any gene code of more than this. 
nPixels = sum(o.HistCounts(:,1,1));
o.z_scoreSHIFT = sum(o.HistValues'.*o.HistCounts)/nPixels;  %Want 0 mean
o.z_scoreSHIFT = zeros(size(o.z_scoreSHIFT));
MeanOfSquare = sum(o.HistValues.^2'.*o.HistCounts)/nPixels;
o.z_scoreSCALE = sqrt(MeanOfSquare - o.z_scoreSHIFT.^2);
SpotColors = (double(o.dpSpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;

%% Get z-scored bleed matrix and bled codes
    
[BleedMatrix,DiagMeasure] = get_bleed_matrix(o,SpotColors,0);
if DiagMeasure<size(o.UseChannels,2)
    error('Bleed matrix not diagonal');
end   
o.z_scoreBleedMatrix = BleedMatrix;
o.z_scoreBledCodes = change_bled_codes(o,o.z_scoreBleedMatrix);
%All BledCodes need same norm for OMP to be fair.
o.z_scoreBledCodes = o.z_scoreBledCodes./vecnorm(o.z_scoreBledCodes(:,:),2,2);

%% Get background codes and append the best ones to get ompBledCodes
[o.BackgroundEigenvectors,o.BackgroundEigenvalues,o.BackgroundMaxGeneDotProduct,...
    o.BackgroundMaxGeneDotProductGene,o.BackgroundEigenvectorTiles] = o.get_background_codes;
o.UseBackgroundEigenvectors = o.BackgroundEigenvalues>=o.BackgroundEigenvalues(...
    min(o.Max_nBackground,length(o.BackgroundEigenvalues))) & ...
    o.BackgroundMaxGeneDotProduct<o.BackgroundMaxGeneDotProductThresh;
o.nBackground = sum(o.UseBackgroundEigenvectors);
o.ompBledCodes = o.z_scoreBledCodes;
nCodes = size(o.z_scoreBledCodes,1);
o.ompBledCodes = zeros(nCodes+o.nBackground,o.nRounds*o.nBP);
o.ompBledCodes(1:nCodes,:) = o.z_scoreBledCodes;
o.ompBledCodes(nCodes+1:nCodes+o.nBackground,:) = o.BackgroundEigenvectors(o.UseBackgroundEigenvectors,:);
end

