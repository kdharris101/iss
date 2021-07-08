%Want to see want determines intensity of spots bleed matrix found from
%i.e. round, channel, spot intensity, gene number etc.
%% First record median,mean, stdev color in each round and channel for good spots 
% found by pixel based method
pxUse = o.pxSpotScore>200;
p = o.BledCodesPercentile;
SpotColors = double(o.pxSpotColors(pxUse,:,:))./p;
m = permute(squeeze(squeeze(SpotColors)),[1 3 2]);
m = squeeze(reshape(m,[],size(m,1)*o.nRounds,o.nBP));
mNorm = bsxfun(@rdivide, m, sqrt(sum(m.^2,2)));
v = eye(o.nBP);
score = mNorm*v';
[~, Channel] = max(score,[],2); % find best cluster for each point
Round = repelem((1:o.nRounds)',length(Channel)/o.nBP,1);
crColor = m(sub2ind(size(m),(1:length(m))',Channel));
%pxData(b,r) = [median,mean,stdev]
pxData = zeros(o.nBP,o.nRounds,3);
for r=1:o.nRounds
    for b=1:o.nBP
        CurrentColor = crColor(Channel==b&Round==r);
        pxData(b,r,1) = median(CurrentColor*p(:,b,r));
        pxData(b,r,2) = mean(CurrentColor*p(:,b,r));
        pxData(b,r,3) = std(CurrentColor*p(:,b,r));
    end
end

%% Get same data for current method of obtaining bleed matrix
%Single BleedMatrix
p = o.BledCodesPercentile;
SpotColors = o.dpSpotColors(o.dpSpotIsolated,:,:)./p;
m = permute(squeeze(squeeze(SpotColors)),[1 3 2]);
m = squeeze(reshape(m,[],size(m,1)*o.nRounds,o.nBP));
m = m(~any(isnan(m),2),:);
[Channel, v, s2] = ScaledKMeans(m, eye(o.nBP));
Round = repelem((1:o.nRounds)',length(Channel)/o.nBP,1);
NotUsed = Channel==0;
Channel(NotUsed)=1;
crColor = m(sub2ind(size(m),(1:length(m))',Channel));
Channel(NotUsed)=0;
%SingleBMData(r,b) = [median,mean,stdev,EigVectorPred];
CurrentMethod.SingleBMData = zeros(o.nBP,o.nRounds,4);
for b=1:o.nBP
    for r=1:o.nRounds 
        CurrentColor = crColor(Channel==b&Round==r);
        CurrentMethod.SingleBMData(b,r,1) = median(CurrentColor*p(:,b,r));
        CurrentMethod.SingleBMData(b,r,2) = mean(CurrentColor*p(:,b,r));
        CurrentMethod.SingleBMData(b,r,3) = std(CurrentColor*p(:,b,r));
        CurrentMethod.SingleBMData(b,r,4) = v(b,b) * sqrt(s2(b)) * p(:,b,r);
    end
end

CurrentMethod.EachRoundBMData = zeros(o.nBP,o.nRounds,4);
%BleedMatrix for each round
for r=o.UseRounds
    m = squeeze(SpotColors(:,:,r)); % data: nCodes by nBases
    m = m(~any(isnan(m),2),:);
    [Channel, v, s2] = ScaledKMeans(m, eye(o.nBP));
    NotUsed = Channel==0;
    Channel(NotUsed)=1;
    crColor = m(sub2ind(size(m),(1:length(m))',Channel));
    Channel(NotUsed)=0;
    for b=1:o.nBP
        CurrentColor = crColor(Channel==b);
        CurrentMethod.EachRoundBMData(b,r,1) = median(CurrentColor*p(:,b,r));
        CurrentMethod.EachRoundBMData(b,r,2) = mean(CurrentColor*p(:,b,r));
        CurrentMethod.EachRoundBMData(b,r,3) = std(CurrentColor*p(:,b,r));
        CurrentMethod.EachRoundBMData(b,r,4) = v(b,b) * sqrt(s2(b)) * p(:,b,r);
    end
end


%% Now see what values you get from various subsets of anchor spots
% I.e. what is used to find the bleed matrices.
load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
SpotColorsRaw = o.dpSpotColors;
SpotColors = SpotColorsRaw./o.BledCodesPercentile;
SpotIsolated = o.dpSpotIsolated;
m = permute(squeeze(squeeze(SpotColors(:,o.UseChannels,o.UseRounds))),[1 3 2]);
m = squeeze(reshape(m,[],size(m,1)*o.nRounds,o.nBP));
mNorm = bsxfun(@rdivide, m, sqrt(sum(m.^2,2)));
v = eye(o.nBP);
score = mNorm*v';

[TopScore, Channel] = max(score,[],2); % find best cluster for each point
Use = TopScore>0;
nSpots = length(SpotIsolated);
SpotNo = repmat((1:nSpots)',o.nRounds,1);
SpotIsolated = double(repmat(SpotIsolated,o.nRounds,1));
Round = repelem((1:o.nRounds)',nSpots,1);
crColor = m(sub2ind(size(m),(1:length(m))',Channel));   %Color in assigned channel (normalised by percentile).
Tile = double(o.dpLocalTile(SpotNo));
LocalYX = o.dpSpotGlobalYX(SpotNo,:)-o.TileOrigin(Tile,:,o.ReferenceRound);
Dist2TileCentre = vecnorm(LocalYX-o.TileCentre,2,2);
GeneNo = o.pSpotCodeNo(SpotNo);
SpotIntensity = prctile(SpotColorsRaw(SpotNo,:)',47.5*100/49.0)';

Dist2NearestSpot = zeros(size(SpotNo));
nTiles = size(o.TileOrigin,1);
for t=1:nTiles
    for r=1:o.nRounds
        for b=1:o.nBP
            Set = Tile==t & Round==r & Channel==b;
            [~,Dist2NearestSpot(Set)]=knnsearch(AllBaseLocalYX{t,b,r},LocalYX(Set,:));
        end
    end
end

%Remove spots not used to get bleed matrix
Data = [SpotNo,crColor,Channel,Round,Tile,Dist2TileCentre,SpotIsolated,GeneNo,...
    TopScore,SpotIntensity,Dist2NearestSpot];
RoundDistNearestSpotThresh = zeros(o.nRounds,1);
Use = int32(TopScore>0.5);
for r=1:o.nRounds
    RoundDistNearestSpotThresh(r) = max(Dist2NearestSpot(Round==r&crColor>0.5));
    Use = Use+int32(Round==r&Dist2NearestSpot<RoundDistNearestSpotThresh(r));
end
Use = Use>=2;
Data = Data(Use,:);
%plotmatrix(Data(:,2:end),Data(:,2:end));

% figure;
% Colors = colormap(lines(o.nBP));
% Colors = repelem(Colors,length(o.UseRounds),1);
% boxplot_x = [];
% boxplot_g = [];
% i=1;
% for b=1:o.nBP
%     for r=1:o.nRounds
%         bxplt = crColor(Use&Channel==b&Round==r);
%         boxplot_x = [boxplot_x;bxplt];
%         boxplot_g = [boxplot_g;i*ones(size(bxplt))];
%         i=i+1;
%     end
% end
% boxplot(boxplot_x,boxplot_g,'Colors',Colors, 'plotstyle', 'compact');
% set(gca,'xtick',(o.nRounds+1)/2:o.nRounds:o.nBP*o.nRounds,'xticklabel',o.bpLabels);
% xlabel('Channels');

%% Get same data for Use spots
%NewMethodData(r,b) = [median,mean,stdev]
NewMethodData = zeros(o.nBP,o.nRounds,3);
for r=1:o.nRounds
    for b=1:o.nBP
        CurrentColor = crColor(Use&Channel==b&Round==r);
        NewMethodData(b,r,1) = median(CurrentColor*p(:,b,r));
        NewMethodData(b,r,2) = mean(CurrentColor*p(:,b,r));
        NewMethodData(b,r,3) = std(CurrentColor*p(:,b,r));
    end
end

%% Plot results, comparing each method of determining bleed matrix to
% the px median.
GroundTruth = pxData(:,:,1);
boxplot_x = [];
boxplot_xNorm = [];
boxplot_g = [];
Methods = {'Current SingleBM','Current MultipleBM','Median of Current Spots','Median of New Spots'};
i=1;
bxplt = GroundTruth-CurrentMethod.SingleBMData(:,:,4);
bxpltNorm = bxplt./squeeze(p);
boxplot_x = [boxplot_x;abs(bxplt(:))];
boxplot_xNorm = [boxplot_xNorm;abs(bxpltNorm(:))];
boxplot_g = [boxplot_g;i*ones(size(bxplt(:)))];
i=i+1;
bxplt = GroundTruth-CurrentMethod.EachRoundBMData(:,:,4);
bxpltNorm = bxplt./squeeze(p);
boxplot_x = [boxplot_x;abs(bxplt(:))];
boxplot_xNorm = [boxplot_xNorm;abs(bxpltNorm(:))];
boxplot_g = [boxplot_g;i*ones(size(bxplt(:)))];
i=i+1;
bxplt = GroundTruth-CurrentMethod.SingleBMData(:,:,1);
bxpltNorm = bxplt./squeeze(p);
boxplot_x = [boxplot_x;abs(bxplt(:))];
boxplot_xNorm = [boxplot_xNorm;abs(bxpltNorm(:))];
boxplot_g = [boxplot_g;i*ones(size(bxplt(:)))];
i=i+1;
bxplt = GroundTruth-NewMethodData(:,:,1);
bxpltNorm = bxplt./squeeze(p);
boxplot_x = [boxplot_x;abs(bxplt(:))];
boxplot_xNorm = [boxplot_xNorm;abs(bxpltNorm(:))];
boxplot_g = [boxplot_g;i*ones(size(bxplt(:)))];
i=i+1;

figure;
subplot(1,3,1);
Colors = colormap(lines(length(Methods)));
boxplot(boxplot_x,boxplot_g,'Colors',Colors, 'plotstyle', 'compact');
set(gca,'xticklabel',[]);
ylabel('Absolute Difference to GroundTruth');
legend(flip(findall(gca,'Tag','Box')), Methods);
subplot(1,3,2);
boxplot(log10(boxplot_x),boxplot_g,'Colors',Colors, 'plotstyle', 'compact');
set(gca,'xticklabel',[]);
ylabel('log(Absolute Difference to GroundTruth)');
subplot(1,3,3);
boxplot(boxplot_xNorm,boxplot_g,'Colors',Colors, 'plotstyle', 'compact');
set(gca,'xticklabel',[]);
ylabel('Normalized Difference to GroundTruth');


