function [PeakLocalYX,PeakSpotColors,PeakCoefs,PeakNeighbourhoodNonZeros,OriginalTile] = ...
    detect_peak_genes_omp(o,coefs,AllSpotColors,AllLocalYX,t)
%% [PeakLocalYX,PeakSpotColors,PeakCoefs,OriginalTile] = ...
%  detect_peak_genes_omp(o,coefs,GoodSpotColors,GoodLocalYX,t)
%
% This finds the local maxima in omp coefficients for each gene
% 
% Input
% o: iss object
% coefs(S,G) is the weighting of gene G for spot S. coefs(S,74:80) are the
% background weightings. Most are zero.
% AllSpotColors(S,b,r) is the intensity for spot S in channel b, round r.
% S should cover all pixel values that don't go off edge of tile in any b,r.
% AllSpotColors has been z-scored. 
% AllLocalYX(S,:) is the corresponding pixel location.
% t is the current tile of interest
%
% Output
% PeakLocalYX{G} contains the YX position of local maxima of gene G.
% PeakSpotColors{G} contains the corresponding spot colors.
% coefs{G} contains corresponding gene and background weightings
% PeakNeighbourhoodNonZeros{G}: In the neighbourhood defined by 
% se1 = strel('disk', o.PixelDetectRadius); around each peak, this
% specifies the number of pixels that have non zero weightings for gene G.
% OriginalTile{G} = t

%% For each gene, find peaks in omp coef images. Keep these as spots going forward
nCodes = length(o.CharCodes);
PeakSpotColors = cell(nCodes,1);
PeakLocalYX = cell(nCodes,1);
PeakCoefs = cell(nCodes,1);
PeakNeighbourhoodNonZeros = cell(nCodes,1);
OriginalTile = cell(nCodes,1);

GeneIm = zeros(max(AllLocalYX));     %Y index is first in zeros
Ind = sub2ind(size(GeneIm),AllLocalYX(:,1),AllLocalYX(:,2));

fprintf('Tile %d: Finding peaks for gene     ',t);
for GeneNo = 1:nCodes    
    g_num = sprintf('%.6f', GeneNo);
    fprintf('\b\b\b\b%s',g_num(1:4));
    
    %Find local maxima in gene image
    GeneIm(Ind) = coefs(:,GeneNo); 
    Small = 1e-6;
    se1 = strel('disk', o.PixelDetectRadius);     %Needs to be bigger than in detect_spots
    Dilate = imdilate(GeneIm, se1);
    MaxPixels = find(GeneIm + Small >= Dilate & GeneIm > 0);
    
    %Get Indices of Good Global Spot Colors / YX
    PeakInd = find(ismember(Ind,MaxPixels));        %As position in Ind = LogProbOverBackGround Index = Good Index
    nPeaks = length(PeakInd);
    %Save information for that gene
    PeakSpotColors{GeneNo} = AllSpotColors(PeakInd,:,:);
    PeakLocalYX{GeneNo} = AllLocalYX(PeakInd,:);
    PeakCoefs{GeneNo} = coefs(PeakInd,:);    
    OriginalTile{GeneNo} = ones(nPeaks,1)*t;
    %Find number of non zero coefs in neighbourhood around each peak
    NonZeroIm = int8(abs(GeneIm)>0);
    NeighNonZeros = imfilter(NonZeroIm,double(se1.Neighborhood));
    ImagePeakInd = sub2ind(size(GeneIm),AllLocalYX(PeakInd,1),AllLocalYX(PeakInd,2)); %Different from MaxPixels
    PeakNeighbourhoodNonZeros{GeneNo} = NeighNonZeros(ImagePeakInd);
    
    
    if o.Graphics==2
        figure; imagesc(GeneIm); colorbar;
        colormap(gca,bluewhitered);
        hold on
        scatter(AllLocalYX(PeakInd,2),AllLocalYX(PeakInd,1),2,'kx');
        hold off
        title(['Gene',num2str(GeneNo),': ',o.GeneNames{GeneNo}]);
    end
    
    
end
fprintf('\n');
end

