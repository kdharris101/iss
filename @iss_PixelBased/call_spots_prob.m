function [o,LookupTable] = call_spots_prob(o,LookupTableInput)
%% [o, LookupTable] = o.call_spots_prob
% calls spots to codes for in-situ sequencing. Run this after find_spots
%
% o: iss object
% LookupTableInput: Can skip making LookupTable if provide one e.g. if
% running but with different o.ScoreScale
% LookupTable: gives the the probabilities that each spot score is explained 
% by each gene. It saves calculating the probabilities explicitly each time.
%
% produces 
% pSpotCodeNo(Spot): gene index for each spot
% pLogProbOverBackground(Spot): log of probability spot can be explained
% by gene relative to probability it can be explained by background.
% pSpotScore(Spot): pLogProbOverBackground of best gene match relative to
% second best gene match at that location.
% pSpotScoreDev(Spot): standard deviation in spot scores across all genes
% at that location.
% pSpotIntensity(Spot): intensity of the spot. Takes into account
% pSpotCodeNo. Calculated by get_spot_intensity.

%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end

if o.ProbMethod==2 && o.ScoreScale~=1
    warning('o.ProbMethod=2 so changing o.ScoreScale to 1');
    o.ScoreScale = 1;
end

%% Assign prob variables that are equal to dot product method
o.pSpotGlobalYX = o.dpSpotGlobalYX;
o.pSpotColors = o.dpSpotColors;
o.pSpotIsolated = o.dpSpotIsolated;
o.pLocalTile = o.dpLocalTile;

%% Make Bleed Matrices
%Only using channels and rounds given by o.UseChannels and o.UseRounds
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end
    
if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

AllChannels = 1:o.nBP;
IgnoreChannels = setdiff(AllChannels,o.UseChannels);
nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);

%Normalise each colour channel by a percentile as to correct for weaker
%colour channels
p = o.BledCodesPercentile;
NormBleedMatrix = o.BleedMatrix;
if o.Graphics
    figure(98043715); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i);
        imagesc(NormBleedMatrix(:,:,i));
        caxis([0 1]);
        title(sprintf('Round %d', o.UseRounds(i)));
        set(gca, 'xtick', 1:o.nBP);
        set(gca, 'XTickLabel', o.bpLabels(AllChannels));
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'yTickLabel', o.bpLabels(AllChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
    drawnow;
end

%Unnormalise Bleed matrices by multiplying rows by percentiles
BleedMatrix = zeros(o.nBP,o.nBP,nRounds);
for r=1:o.nRounds
    for b=1:o.nBP
        BleedMatrix(b,:,r) = p(:,b,r)*NormBleedMatrix(b,:,r);
    end
end

if o.Graphics
    figure(98043764); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i); 
        imagesc(BleedMatrix(:,:,i)); 
        %caxis([0 1]); 
        title(sprintf('Cycle %d', o.UseRounds(i))); 
        set(gca, 'xtick', 1:o.nBP);
        set(gca, 'XTickLabel', o.bpLabels(AllChannels));
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'yTickLabel', o.bpLabels(AllChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
    drawnow;
%     subplot(2,3,6);
%     caxis([0 1]); 
%     axis off
%     colormap hot
% %     colorbar
end

%Set bad colour channels to 0 for all rounds
BleedMatrix(IgnoreChannels,:,:) = 0;
%o.BleedMatrix = NormBleedMatrix;
o.pBleedMatrix = BleedMatrix;
%BleedMatrix = o.pBleedMatrix;
%% Get BledCodes for each gene
% now load in the code book and apply bleeds to it
%codebook_raw = importdata(o.CodeFile);
%CharCode = codebook_raw.textdata(2:end,5);
%GeneName = codebook_raw.textdata(2:end,3);
GeneName = {};
CharCode = {};
fp = fopen(o.CodeFile, 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName=tmp{1};
CharCode=tmp{2};
fclose(fp);

% bit of a hack to get rid of Sst and Npy (assume always in the end)
nCodes = size(CharCode,1) - nnz(cellfun(@(v) strcmp(v(1:2),'SW'), CharCode));

% put them into object o but without the extras
o.CharCodes=CharCode(1:nCodes);
o.GeneNames=GeneName(1:nCodes);

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r=1:o.nRounds
    if r<=o.nRounds-o.nRedundantRounds
        for c=1:nCodes
            [~, NumericalCode(c,r)] = ismember(CharCode{c}(r), o.bpLabels);
        end
    else
        % redundant round - compute codes automatically
        % find pseudobases for this code
        for c=1:nCodes
            PseudoCode = repmat('0',1,o.nRounds-o.nRedundantRounds);
            for p = 1:length(o.RedundantPseudobases)
                PseudoCode(1,ismember(CharCode{c}, o.RedundantPseudobases{p}))=('0'+p);
            end
            % now match them to the redundant codes
            for cc=1:o.nBP
                rrn = r-o.nRounds+o.nRedundantRounds;
                if ~isempty(regexp(PseudoCode, o.RedundantCodes{rrn,cc}, 'once'))
                    NumericalCode(c,r)=cc;
                end
            end
        end
    end
end

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:nRounds
        if any(AllChannels == NumericalCode(i,o.UseRounds(r))) == 0 continue; end
        BledCodes(i,AllChannels+o.nBP*(r-1)) = BleedMatrix(:, find(AllChannels == NumericalCode(i,o.UseRounds(r))), r);
        UnbledCodes(i,AllChannels(find(AllChannels == NumericalCode(i,o.UseRounds(r))))+o.nBP*(r-1)) = 1;
    end
end

o.pBledCodes = BledCodes;

%% Find probability of each spot to each gene
%Comes from P(f) = P_lambda(lambda)P_hist(f-lambda*g) as f = lambda*g+background

%Load histogram data - background prob distribution
%Need to make hist data symmetric and include all data - 0 in middle
%This assumes -NewValMax < min(o.HistValues).
ModHistCounts = o.SmoothHistCounts;     %Smooth extreme values of HistCounts
[NonZeroValIdx,~,~] = ind2sub(size(ModHistCounts),find(ModHistCounts>0.1));
NewValMax = max(max(o.pSpotColors(:)),o.HistValues(max((NonZeroValIdx))));
o.SymmHistValues = -NewValMax:NewValMax;
nBins = length(o.SymmHistValues);
SymmHistCounts = zeros(nBins,o.nBP,o.nRounds);
LastIdx = find(o.HistValues==o.SymmHistValues(end));
if o.SymmHistValues(1)<o.HistValues(1)
    FirstIdx = find(o.SymmHistValues==o.HistValues(1));
    SymmHistCounts(FirstIdx:end,:,:) = ModHistCounts(1:LastIdx,:,:);
else
    FirstIdx = find(o.HistValues==o.SymmHistValues(1));
    SymmHistCounts(:,:,:) = ModHistCounts(FirstIdx:LastIdx,:,:);
end
HistProbs = SymmHistCounts./sum(SymmHistCounts);
o.HistProbs = (HistProbs+o.alpha)./(1+nBins*o.alpha);

%Get Lambda probability distribution for all genes
x = min(o.pSpotColors(:))-1:max(o.pSpotColors(:))-1;    %subsitiution x=lambda*g, -1 due to matlab indexing
o.ZeroIndex = find(x==0);     %need when looking up conv values
LambdaDistAll = zeros(length(x),nCodes,o.nBP,o.nRounds);

fprintf('\nGetting probability distributions for gene   ');
for GeneNo = 1:nCodes
    if GeneNo<10
        fprintf('\b%d',GeneNo);
    else
        fprintf('\b\b%d',GeneNo);
    end
    BledCode = reshape(o.pBledCodes(GeneNo,:),[o.nBP,o.nRounds]);
    numCharCode = str2double(regexp(cell2mat(o.CharCodes(GeneNo)),'\d','match'))+1;
    
    for b=1:o.nBP
        for r=1:o.nRounds
            g = BledCode(b,r);
            if g<0
                error('Predicted bled code is less than 0');
            end
            if o.ProbMethod == 1
                LambdaDist = gampdf(x,o.GammaShape,g/o.GammaShape);
                %Normalise as for small g might not be normalised due to
                %discrete distribution. 
            elseif o.ProbMethod == 2
                if numCharCode(r)==b
                    %for b/r in CharCodes, expect non zero lambda.
                    %g always >0 in this case
                    LambdaDist = raylpdf(x/g,o.RaylConst);
                else
                    %for b/r not in CharCodes, expect approx zero lambda.
                    LambdaDist = (o.ExpConst/2)*exp(-o.ExpConst*abs(x/g));
                end
            end
            LambdaDistAll(:,GeneNo,b,r) = LambdaDist/sum(LambdaDist);
        end
    end
end
o.LambdaDist = LambdaDistAll;

%Store convolution results as look up table
if nargin<2 || isempty(LookupTableInput)    
    LookupTable = nan(length(x),nCodes,o.nBP,o.nRounds);
    fprintf('\nDoing convolutions for Channel  ');
    for b=1:o.nBP
        fprintf('\b%d',b);
        if ismember(b,o.UseChannels)
            for r=1:o.nRounds
                LookupTable(:,:,b,r)=log(conv2(LambdaDistAll(:,:,b,r),o.HistProbs(:,b,r),'same'));
            end
        end
    end
else
    LookupTable = LookupTableInput;
end
fprintf('\n');

if o.ProbMethod == 1
    %Get background Prob using gamma with g=1 i.e. very small prediction of
    %intensity.
    gBackground = 1.0;
    BackgroundGamma = gampdf(x,o.GammaShape,gBackground/o.GammaShape);
    BackgroundGamma = BackgroundGamma/sum(BackgroundGamma);
elseif o.ProbMethod == 2
    %This just shifts o.HistProbs when convolved with it to where x=0.
    BackgroundGamma = dirac(x);
    BackgroundGamma(BackgroundGamma==inf)=1;
end
o.BackgroundProb = zeros(length(x),o.nBP,o.nRounds);
for b=1:o.nBP
    for r=1:o.nRounds
        o.BackgroundProb(:,b,r) = conv(BackgroundGamma,o.HistProbs(:,b,r),'same');
    end
end

%Get log probs for each spot 
LogProbOverBackground = o.get_LogProbOverBackground(o.pSpotColors,LookupTable);
[LogProbOverBackground,SpotCodeNo] = sort(LogProbOverBackground,2,'descend');

o.pLogProbOverBackground = LogProbOverBackground(:,1);
o.pSpotCodeNo = SpotCodeNo(:,1);
o.pSpotScore = LogProbOverBackground(:,1)-LogProbOverBackground(:,2);
%Store deviation in spot scores - can rule out matches based on a low
%deviation.
o.pSpotScoreDev = std(LogProbOverBackground,[],2);
[o.pSpotIntensity,o.pSpotIntensity2] = o.get_spot_intensity(o.pSpotCodeNo,o.pSpotColors);
if nargin<2 || isempty(LookupTableInput)
    save(fullfile(o.OutputDirectory, ['LookupTable',num2str(o.ProbMethod),'.mat']),'LookupTable','-v7.3');
end

end
