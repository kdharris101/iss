function [LogProbOverBackground,LogProbOverBackgroundMatrix] = get_LogProbOverBackground(o,SpotColors,LookupTable)
%% LogProbOverBackground = o.get_LogProbOverBackground(SpotColors,LookupTable);
% This returns the probability that each spot can be explained by each gene
% relative to probability background can explain it alone.
% o: iss object
% SpotColors(s,b,r): color of spot s in channel b, round r.
% LookupTable: stores convolutions of lambda distributions with background
% histograms.
% LogProbOverBackground(s,g) gives the probabilty that spot s can be
% explained by gene o.GeneNames{g)}.
% LogProbOverBackgroundMatrix(i,g) is the LogProbOverBackground in
% channel/round indicated by i for gene g. This does not take account of
% ScoreScale. 
%Variables needed for summing LogProbabilities from lookup table
nCodes = length(o.CharCodes);
nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);
gChannelIndex = int32(repmat(o.UseChannels,1,nRounds));
gRoundIndex = int32(repelem(o.UseRounds,1,nChans));
ChannelIndex = repmat(gChannelIndex,1,nCodes);
RoundIndex = repmat(gRoundIndex,1,nCodes);
GeneIndex = int32(repelem(1:nCodes,1,nRounds*nChans));

nSpots = size(SpotColors,1);
Verbose = nSpots>1000;
SpotColors = int32(SpotColors); %For indexing everything needs to be int32
LogProbOverBackground = zeros(nSpots,nCodes);

LogProbMultiplier = zeros(nRounds*nChans,nCodes);
if nChans == o.nBP && nRounds == o.nRounds
    for g=1:nCodes
        LogProbMultiplier(:,g) = o.UnbledCodes(g,:);
    end
    %How many squares that contribute:
    NormFactor = double(o.nBP*o.nRounds)/double(o.nRounds+o.ScoreScale*(o.nBP*o.nRounds-o.nRounds));
    %Normalise by this to allow valid comparison
    LogProbMultiplier(LogProbMultiplier==0) = o.ScoreScale;
    LogProbMultiplier = LogProbMultiplier*NormFactor;
else
    %If missing some channels, need all rounds/channels to contribute
    %equal.
    if o.ScoreScale~=1
        error('As missing some channels or rounds, you should make o.ScoreScale=1 so all remaining rounds and channels contribute');
    end
    LogProbMultiplier(:) = 1;
end
if Verbose
    fprintf('Percentage of spot probabilities found:       ');
end
for s=1:nSpots
    sSpotColor = SpotColors(s,sub2ind([o.nBP,o.nRounds],gChannelIndex,gRoundIndex));
    SpotIndex = repmat(o.ZeroIndex-1+sSpotColor,1,nCodes); %-1 due to matlab indexing I think
    %Indices = sub2ind(size(LookupTable),SpotIndex,GeneIndex,ChannelIndex,RoundIndex);
    Indices = SpotIndex + (GeneIndex-1)*size(LookupTable,1) +...
         (ChannelIndex-1)*size(LookupTable,1)*size(LookupTable,2)+...
         (RoundIndex-1)*size(LookupTable,1)*size(LookupTable,2)*size(LookupTable,3);
    BackgroundSpotIndex = o.ZeroIndex-1+sSpotColor;
    BackgroundIndices = BackgroundSpotIndex+...
        (gChannelIndex-1)*size(o.BackgroundProb,1)+...
        (gRoundIndex-1)*size(o.BackgroundProb,1)*size(o.BackgroundProb,2);
    sBackgroundLogProb = log(o.BackgroundProb(BackgroundIndices));
    LogProbMatrix = reshape(LookupTable(Indices),[nRounds*nChans,nCodes]);
    LogProbOverBackground(s,:) = nansum((LogProbMatrix-sBackgroundLogProb').*LogProbMultiplier);
    if mod(s,round(nSpots/100))==0 && Verbose
        Percent = sprintf('%.6f', round(s*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end
if Verbose
    fprintf('\n');
end
if nargout>1 && nSpots==1
    LogProbOverBackgroundMatrix = (LogProbMatrix-sBackgroundLogProb');
end
end

