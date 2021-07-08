%% iss: code for processing of in situ sequencing
% Kenneth D. Harris and Xiaoyan Qian
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% o = iss_GroundTruth; % create structure, default parameters
% % change any parameters you want, and set o.FileN
% o = o.extract_and_filter; % create top-hat filtered tiffs for each tile
% o = o.find_spots; % find spot positions in global coordinates
% o = o.call_spots; % allocate each spot to a gene
% o = o.call_cells; % identify cells
%
% this current file iss.m just contains default parameter values
classdef iss_GroundTruth < iss_PixelBased
    
    properties
        %% Gad stuff
        gtGeneNo;     %gtGeneNo(r,b) is the gene number that is present in channel b, 
                      %round r. If it is 0, that means there is no ground truth data
                      %for round r, channel b. 
        gtRounds;               %Rounds that contain ground truth data.
        gtReferenceChannel;     %gtReferenceChannel(r): Channel in ground truth 
                                %round r to register to. 0 if no ground
                                %truth in round r. 
        gtAnchorChannel;        %gtAnchorChannel(r) is anchor channel in ground 
                                %truth round r. 0 if no ground truth in round r.
        
        gtRawLocalYX;
        gtRawIsolated;
        
        gtRegInfo;      %Contains information like Scores, nMatches, Error,
                        %nPcConvergedImg for PCR to o.gtRound.
        gtBigImFiles;   %gtBigImFiles{r,b} is file name of stitched imaged for 
                        %channel b, round r. 
                        
        %gtTruePositiveSet{r,b} is above gtColorTruePositiveThresh(r,b) in channel
        %b of round r. It is also higher in channel b than in
        %any other channels of round r.
        %Found by get_gtTruePositiveSet.m
        gtColorTruePositiveThresh = 'auto';
        %If gtColorTruePositiveThresh is 'auto', will set
        %gtColorTruePositiveThresh(r,b) as
        %gtColorTruePositiveThreshAutoPrctile of o.gtHistCounts(:,b,r)
        gtColorTruePositiveThreshAutoPrctile = 0.99;
        gtColorTruePositiveMinThresh = 250; %gtColorTruePositiveThresh should be above this.
        gtTruePositiveMaxSpots = 5000;  %At most this number of spots in true positive set. 
        gtTruePositiveSet;
        %To be true positive, must be closer to gt local maxima than
        %gtTruePosMaxSep.
        gtTruePosMaxSep = 6;
        %When dealing with false positives, have lower threshold for peaks
        %in gt round/channels than when considering true positives. 
        gtColorFalsePositiveThresh = 'auto';  
        %If gtColorFalsePositiveThresh is 'auto', will set
        %gtColorFalsePositiveThresh(r,b) as
        %gtColorFalsePositiveThreshAutoPrctile of o.gtHistCounts(:,b,r)
        gtColorFalsePositiveThreshAutoPrctile = 0.97;
        gtColorFalsePositiveMinThresh = 100; %gtColorFalsePositiveThresh should be above this.
        %To be false positive, must be further from gt local maxima than
        %gtFalsePosMinSep.
        gtFalsePosMinSep = 10;     
        
        %HistCounts(:,b,r) is the pixel count corresponding to each value
        %in HistValues for channel b, round r
        gtHistCounts;
        
        %pf_gtColor(s,b,r) is the intensity of spot pfSpotGlobalYX(s,:) in channel b,
        %round r for method pf. pf_gtColor(:,:,ImagingRounds)=NaN.
        dp_gtColor;
        p_gtColor;
        px_gtColor;
        omp_gtColor;
        sp_gtColor;
        
        %% Ground truth spot data
        % gtSpotGlobalYX{r,b}(s,:) is the y,x coordinates of
        % spot s found in channel b,round r. 
        gtSpotGlobalYX;
        
        % gtSpotColors{r,b}(s,b',r') contains spot color on each base
        % and imaging round of every spot found in channel b, round r.
        gtSpotColors;
        
        % gt_gtColor{r,b}(s,b',r')is the intensity of spot gtSpotGlobalYX{r,b}(s,:) 
        % in channel b' of round r'. gt_gtColor{r,b}(:,:,ImagingRounds)=NaN.
        gt_gtColor;
        
        %gtLocalTile{r,b}(s) is the tile gtSpotGlobalYX{r,b}(s,:) was found on
        gtLocalTile;
        
        %% Ground truth pixel probabilities
        %This is cell2mat version of gt_SpotGlobalYX. 
        gt_pxSpotGlobalYX;
        
        %This is cell2mat version of gtSpotColors.
        gt_pxSpotColors;
        
        %This is o.gtCodeNo(r,b) for all gt_SpotGlobalYX(r,b).
        gt_pxSpotCodeNo;
        
        %gt_pxSpotIntensity(s) is the modified spot intensity given by
        %get_spot_intensity.m usin
        gt_pxSpotIntensity;
        
        %gt_pxLogProb(s) is sum(ln(Prob(b',r'))/ln(HistProb(SpotColor(b',r'),b',r')))
        %i.e. probability spot can be explained by gene relative to
        %probability it can be explained by background alone.
        %For channel b ground truth round r.
        gt_pxLogProbOverBackground;        
        
        %gt_pxSpotScore(s) is pLogProb-SecondLargest(pLogProb).
        gt_pxSpotScore;
        
        %pxSpotScoreDev(s) is the standard deviation of the log prob of 
        %spot s for all genes. For channel b ground truth round r. 
        gt_pxSpotScoreDev;
        
        %% Ground truth OMP results   
        %These correspond to gt_pxSpotGlobalYX.
        
        %gt_ompScore is the reduction in error caused by introduction of gene in
        %rounds/channels where there is not already a gene. Given by
        %get_omp_score
        gt_ompSpotScore;
        
        %gt_ompCoef(s,g) is the weighting of spot s for gene g. Most are zero.
        %The last nBackground are background codes and are always non zero.
        gt_ompCoefs;
        
        %Hard to get gt_ompNeighbNonZeros as would have to do whole tile
        %again for each one. 
        
        %% True positive / False positive labels
        %If pf = o.CallMethodPrefix(GeneCallingMethod):
        %pf_gtIdentity{r,b}(s) is 1 if pfSpotGlobalYX(s,:) is a true positive
        %for gene gtGeneNo(r,b), 2 if it is a false positive and 0 otherwise.
        %pf_gtFound{r,b}(s) is 1 if gtSpotGlobalYX{r,b}(s,:) was found by the
        %GeneCallingMethod and 2 if not. 0 if s is not in true positive
        %set.           
        p_gtIdentity;
        p_gtFound;        
        px_gtIdentity;
        px_gtFound;
        omp_gtIdentity;
        omp_gtFound;
        sp_gtIdentity;
        sp_gtFound;
    end
end

