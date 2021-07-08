%% iss: code for processing of in situ sequencing
% Kenneth D. Harris and Xiaoyan Qian
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% o = iss_OMP; % create structure, default parameters
% % change any parameters you want, and set o.FileN
% o = o.extract_and_filter; % create top-hat filtered tiffs for each tile
% o = o.find_spots; % find spot positions in global coordinates
% o = o.call_spots; % allocate each spot to a gene
% o = o.call_cells; % identify cells
%
% this current file iss.m just contains default parameter values
classdef iss_OMP < iss_GroundTruth
    
    properties
        %% OMP method variables
        %Spot colors are z scored by subtracting z_scoreSHIFT and dividing by
        %z_scoreSCALE
        z_scoreSHIFT;
        z_scoreSCALE;
        
        % BleedMatrix used to estimate BledCodes in call_spots_omp.
        % Each round and channel is z-scored.
        z_scoreBleedMatrix
        % z_scoreBledCodes(nCodes, nBP*nRounds): code vectors after modeling
        % crosstalk and z-scoring each round/channel.
        z_scoreBledCodes;
        
        %BackgroundEigenvectors are found from spots on tiles in
        %BackgroundEigenvectorTiles;
        BackgroundEigenvectorTiles;
        
        %BackgroundEigenvectors(i,b,r) is the intensity in channel b, round
        %r for the ith eigenvector of the covariance matrix comprised of pixels and
        %rounds of spots with low dot product with any gene z_scoreBledCodes.
        %BackgroundEigenvalue(i) is the corresponding eigenvalue.
        BackgroundEigenvectors;
        BackgroundEigenvalues;
        
        %BackgroundMaxGeneDotProduct(i) is the absolute dot product of
        %BackgroundEigenvectors(i,:) with
        %z_scoreBledCodes(BackgroundMaxGeneDotProductGene(i),:).
        %This is the largest dot product of any gene with BackgroundEigenvectors(i,:)
        %Idea is that large dot product means it is similar to a gene code
        %so reject any larger than BackgroundMaxGeneDotProductThresh.
        BackgroundMaxGeneDotProduct;
        BackgroundMaxGeneDotProductGene;
        BackgroundMaxGeneDotProductThresh = 0.5;
        
        %BackgroundEigenvectors(UseBackgroundEigenvectors,:) are appended
        %to z_scoreBledCodes to form ompBledCodes. Default is to choose
        %first Max_nBackground eigenvectors, but ignoring any with 
        %BackgroundMaxGeneDotProduct>BackgroundMaxGeneDotProductThresh.
        UseBackgroundEigenvectors;
        nBackground;    %size(UseBackgroundEigenvectors).
        Max_nBackground = 7; %there will be at most this number of background codes.
        ompBledCodes;
        
        %OMP stops when reduction in residual drops below
        %ResidualThresh = ResidualThreshParam*SpotSecondLargestIntensity.
        ResidualThreshParam = 0.175;
        %ResidualThresh is clamped between the two values below.
        ResidualThreshMin = 0.0875;
        ResidualThreshMax = 3.0;
        
        %ompMaxGenes is the maximum number of genes that can be assigned to
        %each pixel.
        ompMaxGenes = 6;
        
        %For quality_threshold:
        %Spots must have either ompSpotIntensity>ompIntensityThresh or
        %ompNeighbNonZeros>ompNeighbThresh or ompSpotScore>ompScoreThresh.
        ompIntensityThresh=1.7891e+03;  %Optimized using PyTorch
        ompIntensity2Thresh = 0;
        ompNeighbThresh=19.9863;        %Optimized using PyTorch
        ompNeighbThresh2=2;
        ompScoreThresh = 7.9784         %Optimized using PyTorch
        
        %% OMP method outputs
        
        % ompSpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        % pixel spots are at different locations.
        ompSpotGlobalYX;
        
        % ompSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        ompSpotColors;
        
        %ompSpotIntensity is the modified spot intensity given by
        %get_spot_intensity.m
        ompSpotIntensity;
                      
        %ompScore is the reduction in error caused by introduction of gene in
        %rounds/channels where there is not already a gene. Given by
        %get_omp_score
        ompSpotScore;
                
        %ompSpotCodeNo is the gene found for each spot
        ompSpotCodeNo;
        
        %ompLocalTile(s) is the tile spot s was found on
        ompLocalTile;
        
        %ompCoef(s,g) is the weighting of spot s for gene g. Most are zero.
        %The last nBackground are background codes and are always non zero.
        ompCoefs;
        %ompNeighbNonZeros(s) is the number of pixels in the surrounding
        %region around spot s that have non zero coefficient for the gene
        %specified by ompSpotCodeNo(s). Min value is 1. Region size
        %specified by o.PixelDetectRadius. For default,
        %o.PixelDetectRadius=4, max value of this is 37.
        ompNeighbNonZeros;
        
    end
end

