function iss_view_prob_no_spot(o, FigNo, Norm, LookupTable, GeneNames, SpotNo, Method)
%% iss_view_prob_no_spot(o, FigNo, Norm, LookupTable, GeneNames)
%
% This function lets you view the spot code, gene code,
% and ln(Prob) for all matches for a chosen location.
%
% o: iss object
% FigNo: figure number (default, current figure)
% Norm: normalization option (described below)
% LookupTable: output from call_spots_prob and is saved in
% o.OutputDirectory. It helps do convolutions more quickly
% GenesNames is a cell of gene names that you want to consider e.g.
% [{'Npy'},{'Pvalb'}]. It is case sensitive.
% SpotNo is number of px Spot interested in or is yx location want to look
% at.
%
% Norm = 1: Raw Colors
% Norm = 2: Normalised by o.SpotNormPrctile in each colour channel and round,
% then if o.CallSpotsCodeNorm='WholeCode', normalise so whole code has L2 norm
% of 1 but if o.CallSpotsCodeNorm='Round', normalise so each round has L2 norm of 1.
% Norm = 3: Normalised by percentile for each color channel across all
% rounds


%%
t=0;
if nargin>=2 && nargin<6
    figure(FigNo);
end

if nargin<3 || isempty(Norm)
    Norm = 1;
end

if nargin<5 || isempty(GeneNames)
    GeneNames = o.GeneNames;
end
GeneNumbers = find(ismember(o.GeneNames,GeneNames));
if isempty(GeneNumbers)
    warning('GeneNames not valid, so using all genes')
    GeneNames = o.GeneNames;
    GeneNumbers = find(ismember(o.GeneNames,GeneNames));
end

if nargin<6 || isempty(SpotNo)
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
elseif length(SpotNo)==1
    if nargin<7 || isempty(Method)
        xy = o.pxSpotGlobalYX(SpotNo,[2,1]);
    else
        xy = o.([o.CallMethodPrefix(Method),'SpotGlobalYX'])(SpotNo,[2,1]);
    end
    t = o.([o.CallMethodPrefix(Method),'LocalTile'])(SpotNo);
elseif length(SpotNo)==2
    xy = [SpotNo(2),SpotNo(1)];
end

%% Get in spot color at this position
%Find tile that the point is on and local centered coordinates in reference round
if t==0
    t = o.get_local_tile([xy(2),xy(1)]);
end
LocalYX = [xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound)-o.TileCentre;
SpotColor = zeros(1,o.nBP,o.nRounds);
for r=1:o.nRounds
    for b=1:o.nBP
        
        rbYX = round([LocalYX,1]*o.D(:,:,t,r,b)+o.TileCentre);
        y0 = rbYX(1);
        x0 = rbYX(2);
        if y0>o.TileSz || y0<1 || x0>o.TileSz || x0<1
            warning('Round %d, channel %d on different tile, setting to 0',r,b)
            continue;
        end
        SpotColor(:,b,r) = int32(imread(o.TileFiles{r,t}, b, 'PixelRegion', {[y0 y0], [x0 x0]}))-o.TilePixelValueShift;
    end
end

%% Get matching gene and score values
[LogProbOverBackground,LogProbOverBackgroundMatrix] = get_LogProbOverBackground(o,SpotColor,LookupTable);
[S.LogProbOverBackground,S.CodeNoAll] = sort(LogProbOverBackground,2,'descend');
nCodes = length(o.CharCodes);
S.ProbMatrices = reshape(LogProbOverBackgroundMatrix',[nCodes,o.nBP,o.nRounds]);

S.GeneRank = find(ismember(S.CodeNoAll,GeneNumbers));
S.CodeNoAll = S.CodeNoAll(S.GeneRank);
nCodesToUse = length(S.CodeNoAll);
S.SecondBestLogProbOverBackground = S.LogProbOverBackground(2);
S.LogProbOverBackground = S.LogProbOverBackground(S.GeneRank);
S.CodeIdx = 1;
S.CodeNo = S.CodeNoAll(S.CodeIdx);
S.SpotScore = S.LogProbOverBackground(S.CodeIdx)-S.SecondBestLogProbOverBackground;
S.SpotScoreDev = std(LogProbOverBackground,[],2);
S.SpotIntensity = o.get_spot_intensity(S.CodeNoAll,repmat(SpotColor,[nCodesToUse,1,1]));
S.SpotColor = SpotColor;


%Different Normalisations
S.sqColor = 'g';
if isempty(Norm) || Norm == 1
    cSpotColor = SpotColor;
    S.cBledCodes = o.pBledCodes;
    S.sqColor = 'r';
elseif Norm == 2
    cSpotColor = double(SpotColor)./o.BledCodesPercentile;
    cSpotColor = cSpotColor/sqrt(sum(cSpotColor(:).^2));
    S.cBledCodes = bsxfun(@rdivide, o.BledCodes, sqrt(sum(o.BledCodes.^2,2)));
elseif Norm == 3
    cSpotColor = double(SpotColor);
    NewBleedMatrix = o.pBleedMatrix;
    for b = 1:o.nBP
        bSpotColors = double(o.pSpotColors(:,b,:));
        p = prctile(bSpotColors(:), o.SpotNormPrctile);
        cSpotColor(:,b,:) = cSpotColor(:,b,:)/p;
        NewBleedMatrix(b,:,:) = o.pBleedMatrix(b,:,:)/p;
    end
    S.cBledCodes = change_bled_codes(o,NewBleedMatrix);
end

MeasuredCode = squeeze(cSpotColor);
S.CodeShape = size(MeasuredCode);

%Get square outlining unbled code
gUnbled = reshape(o.UnbledCodes(S.CodeNo,:,:),S.CodeShape);
S.gSquares = zeros(o.nRounds,4);
for r=1:o.nRounds
    try
        S.gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
    end
end

%So don't have to take o object everywhere.
S.Norm = Norm;
S.nBP = o.nBP;
S.bpLabels = o.bpLabels;
S.nRounds = o.nRounds;
S.GeneNames = o.GeneNames;
S.pScoreThresh = o.pScoreThresh;
S.pLogProbThresh = o.pLogProbThresh;
S.pDevThresh = o.pDevThresh;
S.pIntensityThresh = o.pIntensityThresh;
S.UnbledCodes = o.UnbledCodes;
S.MinAllColors = min(o.pSpotColors(:));
S.MaxAllColors = max(o.pSpotColors(:));
S.LambdaDist = o.LambdaDist;
S.SymmHistValues = o.SymmHistValues;
S.HistProbs = o.HistProbs;
S.BackgroundProb = o.BackgroundProb;
S.BackgroundLambdaDist = o.BackgroundLambdaDist;
S.ZeroIndex = o.ZeroIndex;


try
    clf(430476573)
    S.Fig = figure(430476573);
catch
    S.Fig = figure(430476573);
    set(S.Fig,'Position',[352,57,832,748]);
end
S.ax1 = subplot(3,1,1);
S.SpotImage = imagesc(S.ax1, MeasuredCode); colorbar
%caxis([0 max(MeasuredCode(:))]);
if S.Norm~=1
    colormap(gca,bluewhitered);
end
S.ax1.Title.String = 'Spot Code';
S.ax1.YTick = 1:S.nBP;
S.ax1.YTickLabel = S.bpLabels;
S.ax1.YLabel.String = 'Color Channel';
hold on
for r=1:S.nRounds
    rectangle(S.ax1,'Position',S.gSquares(r,:),'EdgeColor',S.sqColor,'LineWidth',2,'LineStyle',':')
end
hold off

S.ax2 = subplot(3,1,2);
S.GeneImage = imagesc(S.ax2, reshape(S.cBledCodes(S.CodeNo,:), S.CodeShape)); colorbar;
S.GeneImageCaxis = [min(S.cBledCodes(:)),max(S.cBledCodes(:))];
caxis(S.GeneImageCaxis);
if S.Norm~=1
    colormap(gca,bluewhitered);
end
%caxis([0 max(cBledCode(:))]);
S.ax2.Title.String = sprintf('Rank %d: Predicted Code for %s, code #%d',S.GeneRank(S.CodeIdx), S.GeneNames{S.CodeNo}, S.CodeNo);
S.ax2.YTick = 1:S.nBP;
S.ax2.YTickLabel = S.bpLabels;
S.ax2.YLabel.String = 'Color Channel';
hold on
for r=1:S.nRounds
    rectangle(S.ax2,'Position',S.gSquares(r,:),'EdgeColor',S.sqColor,'LineWidth',2,'LineStyle',':')
end
hold off

S.ax3 = subplot(3,1,3);
S.ClickPlot = imagesc(S.ax3,squeeze(S.ProbMatrices(S.CodeNo,:,:))); colorbar;
S.ProbImageCaxis = [min([0;S.ProbMatrices(:)]),max([0;S.ProbMatrices(:)])];
caxis(S.ProbImageCaxis);
colormap(gca,bluewhitered);
%caxis([min(ProbMatrix(:)) max(ProbMatrix(:))]);
S.ax3.YTick = 1:S.nBP;
S.ax3.YTickLabel = S.bpLabels;
S.ax3.YLabel.String = 'Color Channel';
S.ax3.XLabel.String = 'Round';
S.ax3.Title.Interpreter = 'Latex';
S.ax3.Title.FontSize = 15;
S.ax3.Title.String = '$\ln\left({\frac{P(spot\,\mid \,gene\,\, and\,\, background)}{P(spot\,\mid \,background)}}\right)$';
%title(sprintf('Log Probability the spot can be explained by gene - Log Probability it can be explained by background alone'));
set(S.ClickPlot,'ButtonDownFcn',{@getCoord,S});
hold on
for r=1:S.nRounds
    rectangle(S.ax3,'Position',S.gSquares(r,:),'EdgeColor','g','LineWidth',2,'LineStyle',':')
end
hold off

S = getFigureTitle(S);
fprintf('yx=(%d,%d): code %d, %s\n', ...
    xy(2),xy(1),...
    S.CodeNo, S.GeneNames{S.CodeNo});

if nCodesToUse>1
    S.sl = uicontrol('style','slide',...
        'unit','pix',...
        'position',[60 20 713 20],...
        'min',1,'max',nCodesToUse,'val',1,...
        'sliderstep',[1/(nCodesToUse-1) 1/(nCodesToUse-1)],...
        'callback',{@sl_call,S});
    set( findall( S.Fig, '-property', 'Units' ), 'Units', 'Normalized' )
end

end

function S = getFigureTitle(S)
%Color different parameters depending if over threshold
if S.SpotScore>S.pScoreThresh
    c1 = [0,0.7,0]; else; c1 = [0,0,0];end
if S.LogProbOverBackground(S.CodeIdx)<S.pLogProbThresh
    c2 = [1,0,0]; else; c2 = [0,0,0];end
if S.SpotScore+S.SpotScoreDev<S.pDevThresh
    c3 = [1,0,0]; else; c3 = [0,0,0];end
if S.SpotIntensity(S.CodeIdx)<S.pIntensityThresh
    c4 = [1,0,0]; else; c4 = [0,0,0];end


S.figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
S.figtitle.String = sprintf('%s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProbOverBackground = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
    '\color[rgb]',c1,S.SpotScore,'\color[rgb]',c2, S.LogProbOverBackground(S.CodeIdx),...
    '\color[rgb]',c3,S.SpotScoreDev,'\color[rgb]',c4,S.SpotIntensity(S.CodeIdx));
%figtitle.Color='red';
drawnow
end


function getCoord(aH,evnt,S)
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
f = S.SpotColor(:,b,r);
x = min(f,S.MinAllColors-1):max(f,S.MaxAllColors-1);

XLim = [min(x)-2000,max(x)+2000];
NormYLim = [round(min(log(S.HistProbs(:)))-4),0];
AltYLim = [0,max([max(S.LambdaDist(:)),max(S.HistProbs(:))])+0.05];

if strcmp(click_type,'normal')
    LogProbPlot = log(conv(S.LambdaDist(:,S.CodeNo,b,r),S.HistProbs(:,b,r),'same'));
    %Get background too
    BackgroundProb = log(S.BackgroundProb(:,b,r));   
    figure(35458);
    plot(x,BackgroundProb,'Color',[0.25, 0.25, 0.25],'LineWidth',0.1);
    
    %set(get(get(P1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on
    plot(x,LogProbPlot,'Color',[0, 0.4470, 0.7410],'LineWidth',0.8);
    xline(f,'-','DisplayName','Spot Value','Color','red','LineWidth',1);   %Or xline(x(o.ZeroIndex-1+f))
    hold off
    leg1 = legend('$P(s\,\mid \,background)$','$P(s\,\mid \,gene\,\, and\,\, background)$','Spot Value');
    set(leg1,'Interpreter','latex');
    %legend('show');
    xlabel('Spot Intensity, $s$','interpreter','latex','FontSize',13)
    ylabel('Log Probability','interpreter','latex','FontSize',13);
    title('Probability distribution for '+string(S.GeneNames(S.CodeNo))+ ', round '+string(r)+' color channel '+string(b-1));
    xlim(XLim);
    ylim(NormYLim);
    
elseif strcmp(click_type,'alt')
    HistZeroIndex = find(S.SymmHistValues == 0);
    x2 = x(x<HistZeroIndex+f);      %So ensure indices>0
    hIndices = HistZeroIndex+f-x2;
    Use = hIndices<length(S.SymmHistValues);
    HistDist = S.HistProbs(hIndices(Use),b,r);
    LambdaIndices = find(x<HistZeroIndex+f);
    figure(9264992);
    plot(x(LambdaIndices(Use)),S.LambdaDist(LambdaIndices(Use),S.CodeNo,b,r),'Color',[0, 0.4470, 0.7410]);
    hold on
    plot(x(LambdaIndices(Use)),HistDist,'Color','red');
    plot(x(LambdaIndices(Use)),S.BackgroundLambdaDist(LambdaIndices(Use)),'Color',[0.25, 0.25, 0.25]);
    hold off
    title({'For spot s, gene g and background distribution $$P_b$$; given $$x=\lambda g$$:',...
        '$$P(s\mid g) = \int P(\lambda)P_b(s-\lambda g)d\lambda = \frac{1}{g}\sum_{x} P\left(\frac{x}{g}\right)P_b(s-x)$$'},...
        'interpreter','latex','FontSize',13)
    xlabel('$x$','interpreter','latex','FontSize',13)
    ylabel('$Probability$','interpreter','latex','FontSize',13);
    leg1 = legend('$$\frac{1}{g}P\left(\frac{x}{g}\right)$$','$P_b(s-x)$','$$\frac{1}{1}P\left(\frac{x}{1}\right)$$');
    set(leg1,'Interpreter','latex');
    set(leg1,'FontSize',11);
    xlim(XLim);
    ylim(AltYLim);      
end

end



function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
S.CodeIdx = round(get(h,'value'));
S.CodeNo = S.CodeNoAll(S.CodeIdx);
S.SpotScore = S.LogProbOverBackground(S.CodeIdx)-S.SecondBestLogProbOverBackground;

S.GeneImage = imagesc(S.ax2, reshape(S.cBledCodes(S.CodeNo,:), S.CodeShape));
colorbar(S.ax2);
caxis(S.ax2,S.GeneImageCaxis);
if S.Norm~=1
    set(S.Fig,'CurrentAxes',S.ax2);
    colormap(S.ax2,bluewhitered);
end
S.ax2.Title.String = sprintf('Rank %d: Predicted Code for %s, code #%d',S.CodeIdx, S.GeneNames{S.CodeNo}, S.CodeNo);
S.ax2.YTick = 1:S.nBP;
S.ax2.YTickLabel = S.bpLabels;
S.ax2.YLabel.String = 'Color Channel';


S.ClickPlot = imagesc(S.ax3,squeeze(S.ProbMatrices(S.CodeNo,:,:)));
colorbar(S.ax3);
caxis(S.ax3,S.ProbImageCaxis);
set(S.Fig,'CurrentAxes',S.ax3);  %Need to set current axis for bluewhitered to work. 
%axes(S.ax3);    
colormap(S.ax3,bluewhitered);
S.ax3.YTick = 1:S.nBP;
S.ax3.YTickLabel = S.bpLabels;
S.ax3.YLabel.String = 'Color Channel';
S.ax3.XLabel.String = 'Round';
S.ax3.Title.Interpreter = 'Latex';
S.ax3.Title.FontSize = 15;
S.ax3.Title.String = '$\ln\left({\frac{P(spot\,\mid \,gene\,\, and\,\, background)}{P(spot\,\mid \,background)}}\right)$';
set(S.ClickPlot,'ButtonDownFcn',{@getCoord,S});

h = findobj('type','rectangle');
delete(h)
%Get square outlining unbled code
gUnbled = reshape(S.UnbledCodes(S.CodeNo,:,:),S.CodeShape);
S.gSquares = zeros(S.nRounds,4);
for r=1:S.nRounds
    try
        S.gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
    end
end
ax_index = 1;
for ax = [S.ax1,S.ax2,S.ax3]
    if ax_index<3
        ax_color = S.sqColor;
    else
        ax_color = 'g';
    end
    hold on
    for r=1:S.nRounds
        rectangle(ax,'Position',S.gSquares(r,:),'EdgeColor',ax_color,'LineWidth',2,'LineStyle',':');
    end
    hold off
    ax_index = ax_index+1;
end
S = getFigureTitle(S);
drawnow

end
