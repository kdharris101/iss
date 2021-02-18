%% Get Spots to use
Use = o.pxSpotScore>220;        %Need good spots but range of spot intensities.
nSpots = sum(Use);
GeneNo = o.pxSpotCodeNo(Use);
SpotColors = double(o.pxSpotColors(Use,:,:));
SpotIntensity = prctile(SpotColors(:,:)',47.5*100/49.0)';

%% Only use Unbled Codes
nCodes = length(o.CharCodes);
LogProbMultiplier = zeros(nCodes,o.nBP,o.nRounds);
for g=1:nCodes
    LogProbMultiplier(g,:,:) = reshape(o.UnbledCodes(g,:),[o.nBP,o.nRounds]);
end
UnbledSpotColors = SpotColors.*LogProbMultiplier(GeneNo,:,:);
BledCodes = reshape(o.pBledCodes,[nCodes,o.nBP,o.nRounds]);
UnbledCodes = BledCodes.*LogProbMultiplier(1:nCodes,:,:);

%% Find predicted codes found using OMP and channels un-normalized
PredScale = zeros(nSpots,2);
fprintf('Percentage of spot coefs found:       ');
for s=1:nSpots
    PredScale(s,1) = omp_specify_atoms(UnbledCodes(GeneNo(s),:)',UnbledSpotColors(s,:)',1);
    if mod(s,round(nSpots/100))==0
        Percent = sprintf('%.6f', round(s*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end

%% Find predicted codes found using OMP and channels normalized
UnbledSpotColors = UnbledSpotColors./o.BledCodesPercentile;
UnbledCodes = UnbledCodes./o.BledCodesPercentile;
fprintf('Percentage of spot coefs found:       ');
for s=1:nSpots
    PredScale(s,2) = omp_specify_atoms(UnbledCodes(GeneNo(s),:)',UnbledSpotColors(s,:)',1);
    if mod(s,round(nSpots/100))==0
        Percent = sprintf('%.6f', round(s*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end
UnbledCodes = BledCodes.*LogProbMultiplier(1:nCodes,:,:);
UnbledSpotColors = SpotColors.*LogProbMultiplier(GeneNo,:,:);
PredScaleUse = 1;       %1 if channels not normalized, 2 if normalized
PredCodes = PredScale(:,PredScaleUse).*UnbledCodes(GeneNo,:,:);
Residual = UnbledSpotColors-PredCodes;
Residual(LogProbMultiplier(GeneNo,:,:)==0)=nan;

%% Illustrate bled code predictions are too weak
figure; 
xline(mean(PredScale(:,PredScaleUse)),'Color','red','LineWidth',1);
hold on
histogram(PredScale(:,PredScaleUse),-0.2:0.1:4,'Normalization','probability');
hold off
xlabel('Scale of BledCode');
ylabel('Fraction of Spots');
legend('Mean');

%Now look at perfect scales in each round and channel
PerfectScale = UnbledSpotColors./UnbledCodes(GeneNo,:,:);
boxplot_x = [];
boxplot_g = [];
i=1;
for b=1:o.nBP
    for r=1:o.nRounds
        IgnoreNan = ~isnan(PerfectScale(:,b,r));
        bxplt = PerfectScale(IgnoreNan,b,r);
        boxplot_x = [boxplot_x;bxplt];
        boxplot_g = [boxplot_g;i*ones(size(bxplt))];
        i=i+1;
    end
end
Colors = colormap(lines(o.nBP));
Colors = repelem(Colors,length(o.UseRounds),1);
figure;
boxplot(boxplot_x,boxplot_g,'Colors',Colors, 'plotstyle', 'compact');
set(gca,'xtick',(o.nRounds+1)/2:o.nRounds:o.nBP*o.nRounds,'xticklabel',o.bpLabels);
xlabel('Channels');
ylabel('Perfect Scale');
yline(1,'LineWidth',1);

%% Analyse influence of high intensity spots on bleed matrix
BleedMatrixOrig = o.BleedMatrix;
SpotIntensity = prctile(o.dpSpotColors(:,:)',49*100/49.0)';     
NormSpotColors = o.dpSpotColors./o.BledCodesPercentile;
%ModSpotColors sets highest value of that spot to 0.
ModSpotColors = NormSpotColors(:,:);
[a,b] = max(ModSpotColors(:,:),[],2);
ModInd = sub2ind(size(ModSpotColors),(1:size(ModSpotColors,1))',b);
ModSpotColors(Ind) = nan;
ModSpotColors = reshape(ModSpotColors,size(ModSpotColors,1),7,7);

%Remove highest intensity value in spots that are likely to be a gene and
%are in the top 10% of intensity values. I.e. see effect on bleed matrix of
%remove large intensity values.
Set = o.dpSpotIsolated&o.dpSpotScore>0.8;
HighIntensityThresh = prctile(SpotIntensity(Set),90);
HighSet = Set & SpotIntensity>HighIntensityThresh;
UseSpotColors = NormSpotColors;
UseSpotColors(HighSet,:,:) = ModSpotColors(HighSet,:,:);
[BleedMatrixHigh,~] = get_bleed_matrix(o,UseSpotColors,o.dpSpotIsolated);
%Remove highest intensity value in spots that are likely to be a gene and
%are in the bottom 10% of intensity values. I.e. see effect on bleed matrix of
%remove small intensity values.
LowIntensityThresh = prctile(SpotIntensity(Set),10);
LowSet = Set & SpotIntensity<LowIntensityThresh;
UseSpotColors = NormSpotColors;
UseSpotColors(LowSet,:,:) = ModSpotColors(LowSet,:,:);
[BleedMatrixLow,~] = get_bleed_matrix(o,UseSpotColors,o.dpSpotIsolated);

figure;
ax1 = subplot(1,2,1);
imagesc(BleedMatrixHigh(:,:,1)-BleedMatrixOrig(:,:,1));
caxis([-0.05,0.05]);
colormap(gca,bluewhitered);
title('BleedMatrix Excluding High Intensities - Original');
ax2 = subplot(1,2,2);
imagesc(BleedMatrixLow(:,:,1)-BleedMatrixOrig(:,:,1));
caxis([-0.05,0.05]);
colormap(gca,bluewhitered);
title('BleedMatrix Excluding Low Intensities - Original');
colorbar;

pBleedMatrixHigh = zeros(o.nBP,o.nBP,o.nRounds);
pBleedMatrixLow = zeros(o.nBP,o.nBP,o.nRounds);
for r=1:o.nRounds
    for b=1:o.nBP
        pBleedMatrixHigh(b,:,r) = o.BledCodesPercentile(:,b,r)*BleedMatrixHigh(b,:,r);
        pBleedMatrixLow(b,:,r) = o.BledCodesPercentile(:,b,r)*BleedMatrixLow(b,:,r);
    end
end
BledCodesHigh = change_bled_codes(o,pBleedMatrixHigh);
BledCodesLow = change_bled_codes(o,pBleedMatrixLow);

g=1;
figure;
ax1 = subplot(1,2,1);
imagesc(reshape(BledCodesHigh(g,:)-o.pBledCodes(g,:),7,7));
caxis([-980,41]);
colormap(gca,bluewhitered);
title('BledCode Excluding High Intensities - Original');
ax2 = subplot(1,2,2);
imagesc(reshape(BledCodesLow(g,:)-o.pBledCodes(g,:),7,7));
caxis([-980,41]);
colormap(gca,bluewhitered);
title('BledCode Excluding Low Intensities - Original');
colorbar;
sgtitle(o.GeneNames{g});


%% Analyse range of residuals for each spot. 
%Expect small range if each spot is scaled version of bled codes

ResidualRange = range(Residual(:,:),2);
[SortPredScale,Index] = sort(PredScale(:,PredScaleUse));
SortResidualRange = ResidualRange(Index);
AvRange = movmedian(SortResidualRange,1000);
figure; scatter(PredScale(:,PredScaleUse),ResidualRange);
hold on
plot(SortPredScale,AvRange,'LineWidth',5);
hold off
xlabel('Scale of BledCode');
ylabel('Residual Range over 7 unbled rounds/channels');

%%
NormChannelResidual = Residual./o.BledCodesPercentile;
boxplot_x = [];
boxplot_g = [];
i=1;
for b=1:o.nBP
    for r=1:o.nRounds
        IgnoreNan = ~isnan(Residual(:,b,r));
        bxplt = Residual(IgnoreNan,b,r);
        boxplot_x = [boxplot_x;bxplt];
        boxplot_g = [boxplot_g;i*ones(size(bxplt))];
        i=i+1;
    end
end
Colors = colormap(lines(o.nBP));
Colors = repelem(Colors,length(o.UseRounds),1);
figure;
boxplot(boxplot_x,boxplot_g,'Colors',Colors, 'plotstyle', 'compact');
set(gca,'xtick',(o.nRounds+1)/2:o.nRounds:o.nBP*o.nRounds,'xticklabel',o.bpLabels);
xlabel('Channels');
    