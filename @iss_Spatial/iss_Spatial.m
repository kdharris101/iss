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
classdef iss_Spatial < iss_OMP
    
    properties       
       %% Spatial method variables
       %spBledCodes are the combined gene and background codes used in
       %call_spots_spatial
       spBledCodes;
       
       %When finding isolated spots to get spot shape from, consider all
       %spots with dot product with any gene bled code above
       %spShapeThresh and further than spShapeIsolationDist from
       %nearest neighbour.
       spShapeThresh = 5;
       spShapeIsolationDist = 10;
       
       %spShapeEigenvectors(:,n) is the spShapeDiamxspShapeDiam size
       %image that best characterizes a spot in the dot product image with
       %eigenvalue spShapeEigenvalues(n).
       %spShapeRad=(spShapeDiam-1)/2
       spShapeDiam;
       spShapeRad;
       spShapeEigenvectors;
       spShapeEigenvalues;
       
       %When finding genes, will take convolution of each gene dot product
       %image with the first nShapeEigUse shape eigenvectors and sum (As
       %orthogonal).
       nShapeEigUse = 5;
       
       %sp_nIter is the number of iterations of convolutions to be carried
       %out on each tile. 
       sp_nIter = 25;
       
        %% Spatial method outputs
        
        % spSpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        % pixel spots are at different locations.
        spSpotGlobalYX;
        
        % spSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        spSpotColors;
        
        %spSpotIntensity is the modified spot intensity given by
        %get_spot_intensity.m
        spSpotIntensity;
        
        %spScore(s) is the biggest residue reduction in shape around
        %spot spSpotGlobalYX(s,:) as a result of removing gene spSpotCodeNo(s).
        spSpotScore;
                
        %spSpotCodeNo is the gene found for each spot
        spSpotCodeNo;
        
        %spLocalTile(s) is the tile spot s was found on
        spLocalTile;
        
        %spCoef(s,a) is the weighting of spot shape spShapeEigenvectors(:,a),
        %gene spSpotCodeNo(s) at location spSpotGlobalYX(s,:).
        spCoefs;
        
        %spSpotIter is the iteration that the particular spot was found in
        %call_spot_spatial;
        spSpotIter;
        
        
    end
end

