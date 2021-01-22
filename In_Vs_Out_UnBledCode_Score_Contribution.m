%% Analysis of particular spot to see contribution to final score 
% from those rounds/channels in Unbled code and those that aren't

SpotNo = 153485;
SpotColor = o.cSpotColors(SpotNo,:,:);

nCodes = length(o.CharCodes);
SpotIndex = repmat(o.ZeroIndex-1+SpotColor(1,:),1,nCodes); %-1 due to matlab indexing I think
gChannelIndex = repmat(1:o.nBP,1,o.nRounds);
gRoundIndex = repelem(1:o.nRounds,1,o.nBP);
ChannelIndex = repmat(gChannelIndex,1,nCodes);
RoundIndex = repmat(gRoundIndex,1,nCodes);
GeneIndex = repelem(1:nCodes,1,o.nRounds*o.nBP);
Indices = sub2ind(size(LookupTable),SpotIndex,GeneIndex,ChannelIndex,RoundIndex);
LogProb_rb = reshape(LookupTable(Indices),[o.nRounds*o.nBP,nCodes]);

HistZeroIndex = find(o.SymmHistValues == 0);
BackgroundIndices = sub2ind(size(o.HistProbs),HistZeroIndex+SpotColor(1,:),gChannelIndex,gRoundIndex);
BackgroundLogProb_rb = log(o.HistProbs(BackgroundIndices));
ProbMatrices = reshape(LogProb_rb',nCodes,o.nBP,o.nRounds)-...
    reshape(BackgroundLogProb_rb,1,o.nBP,o.nRounds);

InCodeMean = zeros(nCodes,1);
OutCodeMean = zeros(nCodes,1);
for g=1:nCodes
    gUnbled = reshape(o.UnbledCodes(g,:,:),[1,7,7]);
    ProbMatrix = ProbMatrices(g,:,:);
    InCodeMean(g) = sum(sum(ProbMatrix.*gUnbled))/7.0;
    OutCodeMean(g) = sum(sum(ProbMatrix.*~gUnbled))/42.0;
end

GadNo = 22;
SncgNo = 63;
yyaxis left
ToPlot = true(nCodes,1);
ToPlot(GadNo) = false;
ToPlot(SncgNo) = false;
scatter(ones(size(InCodeMean(ToPlot))),InCodeMean(ToPlot)*7);
hold on;
scatter(1,InCodeMean(GadNo)*7,200,'x');
scatter(1,InCodeMean(SncgNo)*7,200,'s');
yyaxis right
scatter(ones(size(OutCodeMean(ToPlot)))+1,OutCodeMean(ToPlot)*42);
scatter(2,OutCodeMean(GadNo)*42,200,'x');
scatter(2,OutCodeMean(SncgNo)*42,200,'s');
yyaxis left
scatter(ones(size(OutCodeMean(ToPlot)))+2,InCodeMean(ToPlot)*7+OutCodeMean(ToPlot)*42);
scatter(3,InCodeMean(GadNo)*7+OutCodeMean(GadNo)*42,200,'x');
scatter(3,InCodeMean(SncgNo)*7+OutCodeMean(SncgNo)*42,200,'s');
xticks([1,2,3]);
xlim([0,4]);
xticklabels({'In Unbled Code';'Not in Unbled Code';'Sum'});
title(['How LogProbOverBackground for all genes is distributed for spot ',num2str(SpotNo)]);
