function [PeakLocalYX,PeakSpotColors,PeakLogProbOverBackground,...
    Peak2ndBestLogProb,PeakScoreDev,OriginalTile] = ...
    detect_peak_genes(o,LookupTable,GoodSpotColors,GoodLocalYX,t)
%% [PeakLocalYX,PeakSpotColors,PeakLogProbOverBackground,...
%    Peak2ndBestLogProb,PeakScoreDev,OriginalTile] = ...
%    detect_peak_genes(o,LookupTable,GoodSpotColors,GoodLocalYX,t)
%
% This finds the local maxima in log probability for each gene
% 
% Input
% o: iss object
% LookupTable(s,G,b,r) is the Log probability for spot intensity s-o.ZeroIndex+1 for gene
% G in channel b, round r.
% GoodSpotColors(S,b,r) is the intensity for spot S in channel b, round r.
% S should cover all pixel values that don't go off edge of tile in any b,r.
% GoodLocalYX(S,:) is the corresponding pixel location.
% t is the current tile of interest
%
% Output
% PeakLocalYX{G} contains the YX position of local maxima of gene G.
% PeakSpotColors{G} contains the corresponding spot colors.
% PeakLogProbOverBackground{G} contains the corresponding 
% Log Probability relative to the background.
% Peak2ndBestLogProb{G} contains the log probability relative to the
% background for the second best match at that location.
% PeakScoreDev{G} is the standard deviation of log probability across all
% genes at that location.
% OriginalTile{G} = t

%% Get log probs for each spot 
%Variables needed for summing LogProbabilities from lookup table
nCodes = length(o.CharCodes);
nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);

gChannelIndex = int32(repmat(o.UseChannels,1,nRounds));
gRoundIndex = int32(repelem(o.UseRounds,1,nChans));
ChannelIndex = repmat(gChannelIndex,1,nCodes);
RoundIndex = repmat(gRoundIndex,1,nCodes);
GeneIndex = int32(repelem(1:nCodes,1,nRounds*nChans));
HistZeroIndex = find(o.SymmHistValues == 0); 

nSpots = size(GoodSpotColors,1);
LogProb = zeros(nSpots,nCodes);
BackgroundLogProb = zeros(nSpots,1);

fprintf('Tile %d: Percentage of spot probabilities found:       ',t);
for s=1:nSpots
    sSpotColor = GoodSpotColors(s,sub2ind([o.nBP,o.nRounds],gChannelIndex,gRoundIndex));
    SpotIndex = repmat(o.ZeroIndex-1+sSpotColor,1,nCodes); %-1 due to matlab indexing I think
    %Indices = sub2ind(size(LookupTable),SpotIndex,GeneIndex,ChannelIndex,RoundIndex);
     Indices = SpotIndex + (GeneIndex-1)*size(LookupTable,1) +...
         (ChannelIndex-1)*size(LookupTable,1)*size(LookupTable,2)+...
         (RoundIndex-1)*size(LookupTable,1)*size(LookupTable,2)*size(LookupTable,3);
    LogProb(s,:)=sum(reshape(LookupTable(Indices),[nRounds*nChans,nCodes]));
    BackgroundIndices = HistZeroIndex+sSpotColor+...
        (gChannelIndex-1)*size(o.HistProbs,1)+...
        (gRoundIndex-1)*size(o.HistProbs,1)*size(o.HistProbs,2);
    %BackgroundIndices = sub2ind(size(o.HistProbs),HistZeroIndex+GoodSpotColors(s,:),gChannelIndex,gRoundIndex);
    BackgroundLogProb(s) = sum(log(o.HistProbs(BackgroundIndices)));
    if mod(s,round(nSpots/100))==0
        Percent = sprintf('%.6f', round(s*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end
fprintf('\n');

AllLogProbOverBackground = LogProb-BackgroundLogProb;
clearvars LogProb BackgroundLogProb;

%% For each gene, find peaks in probability images. Keep these as spots going forward
PeakSpotColors = cell(nCodes,1);
PeakLocalYX = cell(nCodes,1);
PeakLogProbOverBackground = cell(nCodes,1);
Peak2ndBestLogProb = cell(nCodes,1);
PeakScoreDev = cell(nCodes,1);
OriginalTile = cell(nCodes,1);

GeneIm = zeros(max(GoodLocalYX));     %Y index is first in zeros
Ind = sub2ind(size(GeneIm),GoodLocalYX(:,1),GoodLocalYX(:,2));

fprintf('Tile %d: Finding peaks for gene     ',t);
for GeneNo = 1:nCodes    
    g_num = sprintf('%.6f', GeneNo);
    fprintf('\b\b\b\b%s',g_num(1:4));
    
    %Find local maxima in gene image
    GeneIm(Ind) = AllLogProbOverBackground(:,GeneNo); 
    Small = 1e-6;
    se1 = strel('disk', o.PixelDetectRadius);     %Needs to be bigger than in detect_spots
    Dilate = imdilate(GeneIm, se1);
    MaxPixels = find(GeneIm + Small >= Dilate);
    
    %Get Indices of Good Global Spot Colors / YX
    PeakInd = find(ismember(Ind,MaxPixels));        %As position in Ind = LogProbOverBackGround Index = Good Index
    nPeaks = length(PeakInd);
    %Save information for that gene
    PeakSpotColors{GeneNo} = GoodSpotColors(PeakInd,:,:);
    PeakLocalYX{GeneNo} = GoodLocalYX(PeakInd,:);
    peakPoverB = AllLogProbOverBackground(PeakInd,:);
    PeakLogProbOverBackground{GeneNo} = peakPoverB(:,GeneNo);    
    PeakScoreDev{GeneNo} = std(peakPoverB,[],2);
    
    %Find 2nd best gene so can give score relative to it
    [~,gInd] = max(peakPoverB,[],2);
    peakPoverB(sub2ind(size(peakPoverB),(1:nPeaks)',gInd))=-inf;
    Peak2ndBestLogProb{GeneNo} = max(peakPoverB,[],2);
    %SortProb = sort(peakPoverB,2,'descend');
    %Peak2ndBestLogProb{GeneNo} = SortProb(PeakInd,2);
    
    
    OriginalTile{GeneNo} = ones(nPeaks,1)*t;
end
fprintf('\n');
end

