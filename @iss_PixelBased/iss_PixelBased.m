%% iss: code for processing of in situ sequencing
% Kenneth D. Harris and Xiaoyan Qian
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% o = iss_PixelBased; % create structure, default parameters
% % change any parameters you want, and set o.FileN
% o = o.extract_and_filter; % create top-hat filtered tiffs for each tile
% o = o.find_spots; % find spot positions in global coordinates
% o = o.call_spots; % allocate each spot to a gene
% o = o.call_cells; % identify cells
%
% this current file iss.m just contains default parameter values
classdef iss_PixelBased < iss_Base
    
    properties
        %% parameters: Probability Method for spot calling
        % BleedMatrix used to estimate BledCodes in call_spots_prob. Unnormalised.
        pBleedMatrix;
        
        %alpha is used for regularisation so don't have any bins with -Inf
        %log probability.
        alpha = 1e-20;
        
        %An exponential curve will be fit where HistCounts falls below this.
        SmoothHistCountLimit = 500;   
        
        %The HistCount is smoothed according to
        %SmoothHistCountLimit*exp(SmoothHistCountParam*x);
        %SmoothHistCountParamGuess is the start point to find SmoothHistCountParam.
        SmoothHistCountParamGuess = -0.003;
        
        %SymmHistValues is -max(HistValues(HistCount>0)):max(HistValues(HistCount>0)).
        %Required as needs to be symmetric for the convolution
        SymmHistValues;
        
        %HistProbs(:,b,r) is the probability corresponding to each value
        %in SymmHistValues for channel b, round r.
        %(HistCounts/nPixels+o.alpha)./(1+nBins*o.alpha);
        HistProbs;
        
        % pBledCodes(nCodes, nBP*nRounds): code vectors after modeling
        % crosstalk and un-nomalising.
        pBledCodes; 
        
        %ProbMethod determines Lambda dist, either 1 or 2. 
        ProbMethod = 1;
        
        %LambdaDist(:,g,b,r) is the probability distribution of lambda for
        %gene g, channel b and round r. 
        %ProbMethod = 1: Gamma in all rounds/channels using GammaShape.
        %ProbMethod = 2: Rayleigh if appear in CharCodes using RaylConst,
        %Exp otherwise using ExpConst.
        LambdaDist;           
        
        %BackgroundLambdaDist(:) is the background distribution of lambda
        %for all genes, channels and rounds. It is just the LambdaDist
        %where the predicted bled code is 1 (i.e. a very small value).
        BackgroundLambdaDist;
        
        %GammaShape is the shape parameter used in the gamma distribution for
        %the estimated distribution of lambda such that cSpotColors =
        %Lambda*pBledCode
        GammaShape = 3.0;
        
        %RaylConst is the constant used in the rayleigh distribution for
        %the estimated distribution of lambda such that cSpotColors =
        %Lambda*pBledCode
        %1 means distribution peaks at bled code number. 
        RaylConst = 1.0;
        
        %ExpConst is same as above but exponenital distribution used for
        %all rounds/channels that don't appear in CharCode for each gene.
        ExpConst = 3.5;
        
        %ZeroIndex is the index of the 0 value in
        %min(cSpotColors(:)):max(cSpotColors(:)) used for convolutions.
        %Needed to find values in lookup table.
        ZeroIndex;        
        
        %BackgroundProb(:,b,r) is the reference probability distribution
        %for visualisation - 
        %ProbMethod=1: probability distribution if the bled code had a value of 1.
        %ProbMethod=2: HistProbs(:,b,r).
        BackgroundProb;
        
        %ScoreScale is the contribution each round/channel not in Unbled
        %code contributes to LogProbOverBackground compared to each
        %round/channel in Unbled code.
        %0 means only use 7 rounds/channels in unbled code.
        %1 means use all 49 rounds/channels equally
        %In between uses all rounds/channels but unbled code contributes
        %more.
        %If ProbMethod=2, this is 1. 
        ScoreScale = 0;
        
        %pIntensityThresh is the value pSpotIntensity(s) needs to exceed for spot s
        %to count
        pIntensityThresh = 100;
        pIntensityThresh2 = 50;
        
        %pSpotIntensity2 needs to be greater than pIntensity2Thresh
        pIntensity2Thresh = 0;
        
        %pLogProbThresh is the value pLogProbOverBackground(s) needs to exceed for spot s
        %to count
        pLogProbThresh = 0;
        pLogProbThresh2 = 0;
        
        %pScoreThresh is the value pSpotScore(s) needs to exceed for spot s
        %to count
        pScoreThresh = 80;       
        
        %If pSpotScore(s) < pScoreThresh but pSpotScore(s) > pScoreThresh2
        %and has high intensity then will count as spot.
        pScoreThresh2 = 0;
        
        %A spot must have pSpotScore(s)+pSpotScoreDev(s) > pDevThresh to
        %count - avoid spots with similar score to all genes.
        pDevThresh = 6;
        
        %Values used in quality_threshold for prob method
        pQualThresh1 = -220; %Optimized using PyTorch
        pQualParam1 = 1.5;    %Optimized using PyTorch
        pQualThresh2 = 15.634; %Optimized using PyTorch
        pQualParam2 = 0;    %Optimized using PyTorch
        pQualThresh3 = 0;
        pQualThresh4 = 0;
        
        %% variables: spot calling outputs - prob method
        % pSpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        % pixel spots are at different locations.
        pSpotGlobalYX;      %Same as dpSpotGlobalYX
        
        % pSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        pSpotColors;        %Same as dpSpotColors
        
        % pSpotIsolated(Spot) is a binary array saying if the spot is well isolated
        % again for combinatorial spots only
        pSpotIsolated;      %Same as dpSpotIsolated
        
        % pLocalTile(s) is the tile spot cSpotColors(s,:,:) was found on
        pLocalTile;        %Same as dpLocalTile
        
        %pSpotIntensity is the modified spot intensity given by
        %get_spot_intensity.m
        pSpotIntensity;
        
        %pSpotIntensity2 is the median intensity in the 7 rounds/channels
        %specified by the gene assigned i.e. pSpotCodeNo.
        pSpotIntensity2;
        
        %pLogProb is sum(ln(Prob(b,r))/ln(HistProb(SpotColor(b,r),b,r)))
        %i.e. probability spot can be explained by gene relative to
        %probability it can be explained by background alone.
        pLogProbOverBackground;        
        
        %pSpotScore is for the sorted array pLogProb,
        %pLogProb(1)-pLogProb(2)
        pSpotScore;        
        
        %pSpotScoreDev(s) is the standard deviation of the log prob of spot s
        %for all genes
        pSpotScoreDev;
        
        %pSpotCodeNo is the gene found for each spot
        pSpotCodeNo;
        
         %% PixelBased method parameters       
        
        %PixelFileMaxTiles is approximately the maximum number of tiles
        %that can be stored in a single file. Output data to files so don't
        %get memory problems.
        PixelFileMaxTiles = 6;
        
        %PixelFileNames contains the names of files in which pixel method data
        %is stored
        PixelFileNames;
        
        %PixelDetectRadius is the radius of the filter used to find local
        %maxima in gene images
        PixelDetectRadius = 4;
        
        %Have to do initial filtering so don't have too many points, only
        %keep points with pxSpotScore>pxInitialScoreThresh or
        %pxLogProbOverBackground>pxInitialProbThresh
        pxInitialScoreThresh = 0;   %I.e. only want spots that are at least 2nd best at their position
        pxInitialProbThresh = -5;
        
        %% PixelBased method outputs
        
        % pxSpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        % pixel spots are at different locations.
        pxSpotGlobalYX;
        
        % pxSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        pxSpotColors;
        
        %pxSpotIntensity is the modified spot intensity given by
        %get_spot_intensity.m
        pxSpotIntensity;
        
        %pxLogProb is sum(ln(Prob(b,r))/ln(HistProb(SpotColor(b,r),b,r)))
        %i.e. probability spot can be explained by gene relative to
        %probability it can be explained by background alone.
        pxLogProbOverBackground;        
        
        %pxSpotScore is pLogProb -max(pLogProb(SpotCodeNo~=pSpotCodeNo))
        pxSpotScore;
        
        %pxSpotScoreDev(s) is the standard deviation of the log prob of spot s
        %for all genes
        pxSpotScoreDev;
        
        %pxSpotCodeNo is the gene found for each spot
        pxSpotCodeNo;
        
        %pxLocalTile(s) is the tile spot s was found on
        pxLocalTile;
        
        
    end
end

