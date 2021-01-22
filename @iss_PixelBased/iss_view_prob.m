function SpotNo = iss_view_prob(o, FigNo, Norm, Method, IncludeGT, SpotNum)
%% SpotNo = iss_view_prob(o, FigNo, Norm, Method, SpotNum)
%
% This function lets you view the spot code, gene code
% and ln(Prob) for the best match for a chosen spot.
%
% o: iss object
% FigNo: figure number (default, current figure)
% Norm: normalization option (described below)
% Method: 'Prob' or 'Pixel' to consider gene assignments given
% by o.pSpotCodeNo or o.pxSpotCodeNo respectively.
% IncludeGT: if true, will also plot the ground truth rounds.
% SpotNum: number of spot to analyze (default is to click)
% SpotNo: returns the number of the spot analyzed.
%
% Norm = 1: Raw Colors
% Norm = 2: Normalised by o.SpotNormPrctile in each colour channel and round,
% then if o.CallSpotsCodeNorm='WholeCode', normalise so whole code has L2 norm
% of 1 but if o.CallSpotsCodeNorm='Round', normalise so each round has L2 norm of 1.
% Norm = 3: Normalised by percentile for each color channel across all
% rounds


%%
if nargin>=6
    SpotNo = SpotNum;
else
    if nargin>=2
        figure(FigNo);
    end
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    S = evalin('base', 'issPlot2DObject');
    if nargin<4 || isempty(Method)
        if strcmpi('DotProduct',S.CallMethod)
            Method = 'Prob';
        else
            Method = S.CallMethod;
        end
    elseif ~strcmpi(S.CallMethod,Method)
        S.SpotYX = o.([o.CallMethodPrefix(Method),'SpotGlobalYX']);
        S.QualOK = 1;
    end
    InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
    PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
    [~,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);      
end

if nargin<4 || isempty(IncludeGT)
    IncludeGT = false;
end

if ~strcmpi('Prob',Method) && ~strcmpi('Pixel',Method)
    error('Spot calling method not valid, should be Prob or Pixel');
end
pf = o.CallMethodPrefix(Method);        %Method prefix
%Different parameters for different methods
CodeNo = o.([pf,'SpotCodeNo'])(SpotNo);
SpotColor = o.([pf,'SpotColors'])(SpotNo,:,:);
SpotScore = o.([pf,'SpotScore'])(SpotNo);
LogProbOverBackground = o.([pf,'LogProbOverBackground'])(SpotNo);
SpotScoreDev = o.([pf,'SpotScoreDev'])(SpotNo);
SpotIntensity = o.([pf,'SpotIntensity'])(SpotNo);
SpotGlobalYX = o.([pf,'SpotGlobalYX'])(SpotNo,:);

if nargin<3 || isempty(Norm)
    Norm = 1;
end

%Different Normalisations
if isempty(Norm) || Norm == 1
    cSpotColor = double(SpotColor);
    cBledCodes = o.pBledCodes;
elseif Norm == 2
    cSpotColor = double(SpotColor)./o.BledCodesPercentile;
    Norm2SpotColorScale = sqrt(sum(cSpotColor(:).^2));
    cSpotColor = cSpotColor/Norm2SpotColorScale;
    cBledCodes = bsxfun(@rdivide, o.BledCodes, sqrt(sum(o.BledCodes.^2,2)));
elseif Norm == 3
    cSpotColor = double(SpotColor);
    NewBleedMatrix = o.pBleedMatrix;
    for b = 1:o.nBP
        bSpotColors = double(o.pSpotColors(:,b,:));
        p = prctile(bSpotColors(:), o.SpotNormPrctile);
        cSpotColor(:,b,:) = cSpotColor(:,b,:)/p;
        NewBleedMatrix(b,:,:) = o.pBleedMatrix(b,:,:)/p;
    end
    cBledCodes = change_bled_codes(o,NewBleedMatrix);
end

MeasuredCode = squeeze(cSpotColor);
CodeShape = size(MeasuredCode);
BledCode = cBledCodes(CodeNo,:);
ProbMatrix = get_prob_matrix(o,squeeze(SpotColor),CodeNo);
caxis_lims = [0, max(MeasuredCode(:))];

if IncludeGT
    gtSpotColor = o.([pf,'_gtColor'])(SpotNo,:,o.gtRounds);
    if Norm==3
        for b = 1:o.nBP
            bSpotColors = o.pSpotColors(:,b,:);
            p = double(prctile(bSpotColors(:), o.SpotNormPrctile));
            gtSpotColor(:,b,:) = gtSpotColor(:,b,:)/p;
        end
    elseif Norm==2
        gtSpotColor = double(gtSpotColor)./o.BledCodesPercentile(:,:,1:length(o.gtRounds));
        gtSpotColor = gtSpotColor/Norm2SpotColorScale;
    end
    MeasuredCode = [MeasuredCode,squeeze(gtSpotColor)];
end

%Get square outlining unbled code
gUnbled = reshape(o.UnbledCodes(CodeNo(1),:,:),CodeShape);
gSquares = zeros(o.nRounds,4);
for r=1:o.nRounds
    try
        gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
    end
end

try
    clf(430476533)
    figure(430476533)
catch
    figure(430476533)
end
subplot(3,1,1);
imagesc(MeasuredCode); colorbar
caxis(caxis_lims);
title(sprintf('Spot Code'));
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel');
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor','r','LineWidth',2,'LineStyle',':')
end
if IncludeGT
    set(gca, 'xtick', [1:o.nRounds,o.gtRounds]);
    XLabels = string(1:o.nRounds);
    for i=1:length(o.gtRounds)
        XLabels = [XLabels,["gt"+string(i)]];
    end
    set(gca, 'XTickLabel', XLabels);
    for r=o.gtRounds
        for b=1:o.nBP
            if o.gtGeneNo(r,b)==0; continue; end
            if o.gtGeneNo(r,b)==CodeNo
                rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor','r','LineWidth',2.5,'LineStyle','-');
            else
                rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor','k','LineWidth',2.5,'LineStyle','-');
            end
        end
    end
end
hold off

subplot(3,1,2)
imagesc(reshape(BledCode, CodeShape)); colorbar
%caxis([0 max(cBledCode(:))]);
title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel');
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor','r','LineWidth',2,'LineStyle',':')
end
hold off

ClickPlot = subplot(3,1,3);
ClickPlot(1) = imagesc(ProbMatrix); colorbar
colormap(gca,bluewhitered);
%caxis([min(ProbMatrix(:)) max(ProbMatrix(:))]);
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel');
xlabel('Round');
title('$\ln\left({\frac{P(spot\,\mid \,gene\,\, and\,\, background)}{P(spot\,\mid \,background)}}\right)$','interpreter','latex','FontSize',15)
%title(sprintf('Log Probability the spot can be explained by gene - Log Probability it can be explained by background alone'));
set(ClickPlot,'ButtonDownFcn',{@getCoord,o,SpotNo,CodeNo,SpotColor});
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor','g','LineWidth',2,'LineStyle',':')
end
hold off

%Color different parameters depending if over threshold
if SpotScore>o.pScoreThresh
    c1 = [0,0.7,0]; else; c1 = [0,0,0];end
if LogProbOverBackground<o.pLogProbThresh
    c2 = [1,0,0]; else; c2 = [0,0,0];end
if SpotScore+SpotScoreDev<o.pDevThresh
    c3 = [1,0,0]; else; c3 = [0,0,0];end
if SpotIntensity<o.pIntensityThresh
    c4 = [1,0,0]; else; c4 = [0,0,0];end

set(gcf,'Position',[350 100 1000 850])
figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
figtitle.String = sprintf('%s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProbOverBackground = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
    '\color[rgb]',c1,SpotScore,'\color[rgb]',c2, LogProbOverBackground,...
    '\color[rgb]',c3,SpotScoreDev,'\color[rgb]',c4,SpotIntensity);
%figtitle.Color='red';
drawnow

fprintf('Spot %d at yx=(%d,%d): code %d, %s\n', ...
    SpotNo, SpotGlobalYX(1),SpotGlobalYX(2),...
    CodeNo, o.GeneNames{CodeNo});
end

function getCoord(aH,evnt,o,SpotNo,CodeNo,SpotColor)
%This plots a graph showing the variation of probability with spot
%intensity when a left click is applied on a square in the LogProb plot.
%When a right click is applied, a plot showing the individual distributions
%contributing to the LogProb(r,b) in the LogProb plot appears.
drawnow
fig = ancestor(aH,'figure');
click_type = get(fig,'SelectionType');
ClickLoc = evnt.IntersectionPoint(1:2);
r = round(ClickLoc(1));
b = round(ClickLoc(2));
f = SpotColor(:,b,r);
x = min(f,min(o.pSpotColors(:))-1):max(f,max(o.pSpotColors(:))-1);
XLim = [min(x)-2000,max(x)+2000];
NormYLim = [round(min(log(o.HistProbs(:)))-4),0];
AltYLim = [0,max([max(o.LambdaDist(:)),max(o.HistProbs(:))])+0.05];

if strcmp(click_type,'normal')
    LogProbPlot = log(conv(o.LambdaDist(:,CodeNo,b,r),o.HistProbs(:,b,r),'same'));
    %Get background too
    BackgroundProb = log(o.BackgroundProb(:,b,r));   
    figure(35428);
    plot(x,BackgroundProb,'Color',[0.25, 0.25, 0.25],'LineWidth',0.1);
    
    %set(get(get(P1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on
    plot(x,LogProbPlot,'Color',[0, 0.4470, 0.7410],'LineWidth',0.8);
    xline(f,'-','DisplayName','Spot '+string(SpotNo)+ ' Value','Color','red','LineWidth',1);   %Or xline(x(o.ZeroIndex-1+f))
    hold off
    leg1 = legend('$P(s\,\mid \,background)$','$P(s\,\mid \,gene\,\, and\,\, background)$','Spot '+string(SpotNo)+ ' Value');
    set(leg1,'Interpreter','latex');
    %legend('show');
    xlabel('Spot Intensity, $s$','interpreter','latex','FontSize',13)
    ylabel('Log Probability','interpreter','latex','FontSize',13);
    title('Probability distribution for '+string(o.GeneNames(CodeNo))+ ', round '+string(r)+' color channel '+string(b-1));
    xlim(XLim);
    ylim(NormYLim);
    
elseif strcmp(click_type,'alt')
    HistZeroIndex = find(o.SymmHistValues == 0);
    x2 = x(x<HistZeroIndex+f);      %So ensure indices>0
    hIndices = HistZeroIndex+f-x2;
    Use = hIndices<length(o.SymmHistValues);
    HistDist = o.HistProbs(hIndices(Use),b,r);
    LambdaIndices = find(x<HistZeroIndex+f);
    figure(9264892);
    plot(x(LambdaIndices(Use)),o.LambdaDist(LambdaIndices(Use),CodeNo,b,r));
    hold on
    plot(x(LambdaIndices(Use)),HistDist,'Color','red');
    hold off
    title({'For spot s, gene g and background distribution $$P_b$$; given $$x=\lambda g$$:',...
        '$$P(s\mid g) = \int P(\lambda)P_b(s-\lambda g)d\lambda = \frac{1}{g}\sum_{x} P\left(\frac{x}{g}\right)P_b(s-x)$$'},...
        'interpreter','latex','FontSize',13)
    xlabel('$x$','interpreter','latex','FontSize',13)
    ylabel('$Probability$','interpreter','latex','FontSize',13);
    leg1 = legend('$$\frac{1}{g}P\left(\frac{x}{g}\right)$$','$P_b(s-x)$');
    set(leg1,'Interpreter','latex');
    set(leg1,'FontSize',11);
    xlim(XLim);
    ylim(AltYLim);   
end

end
