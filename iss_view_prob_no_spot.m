function iss_view_prob_no_spot(o, FigNo, Norm, LookupTable, GeneNames)
    %% This function lets you view the spot code, gene code,
    %spotcode-lambda*bledcode and ln(Prob) for all matches for a chosen
    %location.
    %Now have option for different normalisations
    %Norm = 1: Raw Colors
    %Norm = 2: Normalised by o.SpotNormPrctile in each colour channel and round,
    %then if o.CallSpotsCodeNorm='WholeCode', normalise so whole code has L2 norm
    %of 1 but if o.CallSpotsCodeNorm='Round', normalise so each round has L2 norm of 1.
    %Norm = 3: Normalised by percentile for each color channel across all
    %rounds
    %LookupTable is an output from call_spots_prob and is saved in
    %o.OutputDirectory. It helps do convolutions more quickly
    %GenesNames is a cell of gene names that you want to consider e.g.
    %[{'Npy'},{'Pvalb'}]. It is case sensitive.
    
    if nargin>=2
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

    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    
    %% Get in spot color at this position
    %Find tile that the point is on and local centered coordinates in reference round
    Dist2Tiles = [xy(2),xy(1)]-o.TileOrigin(:,:,o.ReferenceRound);
    Dist2Tile = min(sum(Dist2Tiles(Dist2Tiles(:,1)>=0 & Dist2Tiles(:,2)>=0,:),2));
    t = find(sum(Dist2Tiles,2)==Dist2Tile);
    LocalYX = [xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound)-o.TileCentre;
    
    SpotColor = zeros(1,o.nBP,o.nRounds);
    for r=1:o.nRounds       
        for b=1:o.nBP

            rbYX = round(o.A(b)*([LocalYX,1]*o.D(:,:,t,r))+o.TileCentre);
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
    nCodes = length(o.CharCodes);
    gChannelIndex = repmat(1:o.nBP,1,o.nRounds);
    gRoundIndex = repelem(1:o.nRounds,1,o.nBP);
    ChannelIndex = repmat(gChannelIndex,1,nCodes);
    RoundIndex = repmat(gRoundIndex,1,nCodes);
    GeneIndex = repelem(1:nCodes,1,o.nRounds*o.nBP);
    HistZeroIndex = find(o.SymmHistValues == 0);
    
    SpotIndex = repmat(o.ZeroIndex-1+SpotColor(1,:),1,nCodes); %-1 due to matlab indexing I think
    Indices = sub2ind(size(LookupTable),SpotIndex,GeneIndex,ChannelIndex,RoundIndex);
    LogProb_rb = reshape(LookupTable(Indices),[o.nRounds*o.nBP,nCodes]);
    LogProbAll = sum(LogProb_rb);
    BackgroundIndices = sub2ind(size(o.HistProbs),HistZeroIndex+SpotColor(1,:),gChannelIndex,gRoundIndex);
    BackgroundLogProb_rb = log(o.HistProbs(BackgroundIndices));
    BackgroundLogProb = sum(BackgroundLogProb_rb);
    S.ProbMatrices = reshape(LogProb_rb',nCodes,o.nBP,o.nRounds)-...
        reshape(BackgroundLogProb_rb,1,o.nBP,o.nRounds);
    
    [LogProbAll,S.CodeNoAll] = sort(LogProbAll,2,'descend');
    S.GeneRank = find(ismember(S.CodeNoAll,GeneNumbers));
    S.CodeNoAll = S.CodeNoAll(S.GeneRank);
    nCodesToUse = length(S.CodeNoAll);
    S.LogProbOverBackground = LogProbAll-BackgroundLogProb;
    S.SecondBestLogProbOverBackground = S.LogProbOverBackground(2);
    S.LogProbOverBackground = S.LogProbOverBackground(S.GeneRank);
    S.CodeIdx = 1;
	S.CodeNo = S.CodeNoAll(S.CodeIdx);
	S.SpotScore = S.LogProbOverBackground(S.CodeIdx)-S.SecondBestLogProbOverBackground;
	S.SpotScoreDev = std(LogProbAll,[],2);
	S.SpotIntensity = o.get_spot_intensity(S.CodeNoAll,repmat(SpotColor,[nCodesToUse,1,1]));
    S.SpotColor = SpotColor;
    
    
    
    %Different Normalisations
    if isempty(Norm) || Norm == 1
        cSpotColor = SpotColor;
        S.cBledCodes = o.pBledCodes;
    elseif Norm == 2
        if strcmpi(o.BleedMatrixType,'Separate')
            p = prctile(o.cSpotColors, o.SpotNormPrctile);
        elseif strcmpi(o.BleedMatrixType,'Single')
            p = zeros(1,o.nBP,o.nRounds);
            for b = 1:o.nBP
                bSpotColors = o.cSpotColors(:,b,:);
                p(:,b,:) = prctile(bSpotColors(:), o.SpotNormPrctile);
            end
        end
        cSpotColor = SpotColor./p;
        cSpotColor = cSpotColor/sqrt(sum(cSpotColor(:).^2));
        %cSpotColor = o.cNormSpotColors(SpotNo,:,:);
        S.cBledCodes = bsxfun(@rdivide, o.BledCodes, sqrt(sum(o.BledCodes.^2,2)));
        %cBledCodes = o.NormBledCodes;
    elseif Norm == 3
        cSpotColor = SpotColor;
        NewBleedMatrix = o.pBleedMatrix;
        for b = 1:o.nBP
            bSpotColors = o.cSpotColors(:,b,:);
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
    S.nBP = o.nBP;
    S.bpLabels = o.bpLabels;
    S.nRounds = o.nRounds;
    S.GeneNames = o.GeneNames;
    S.pScoreThresh = o.pScoreThresh;
    S.pLogProbThresh = o.pLogProbThresh;
    S.pDevThresh = o.pDevThresh;
    S.pIntensityThresh = o.pIntensityThresh;
    S.UnbledCodes = o.UnbledCodes;
    S.MinAllColors = min(o.cSpotColors(:));
    S.MaxAllColors = max(o.cSpotColors(:));
    S.LambdaDist = o.LambdaDist;
    S.SymmHistValues = o.SymmHistValues;
    S.HistProbs = o.HistProbs;
    
    
    try
        clf(430476573)
        S.Fig = figure(430476573);
    catch
        S.Fig = figure(430476573);
    end
    set(S.Fig,'Position',[350 100 1000 850]);
    S.ax1 = subplot(3,1,1);
    S.SpotImage = imagesc(S.ax1, MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    S.ax1.Title.String = 'Spot Code';
    S.ax1.YTick = 1:S.nBP;
    S.ax1.YTickLabel = S.bpLabels;
    S.ax1.YLabel.String = 'Color Channel';
    hold on
    for r=1:S.nRounds
        rectangle(S.ax1,'Position',S.gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
    end
    hold off
    
    S.ax2 = subplot(3,1,2);
    S.GeneImage = imagesc(S.ax2, reshape(S.cBledCodes(S.CodeNo,:), S.CodeShape)); colorbar;
    %caxis([0 max(cBledCode(:))]);
    S.ax2.Title.String = sprintf('Rank %d: Predicted Code for %s, code #%d',S.GeneRank(S.CodeIdx), S.GeneNames{S.CodeNo}, S.CodeNo);
    S.ax2.YTick = 1:S.nBP;
    S.ax2.YTickLabel = S.bpLabels;
    S.ax2.YLabel.String = 'Color Channel';
    hold on
    for r=1:S.nRounds
        rectangle(S.ax2,'Position',S.gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
    end
    hold off
    
    S.ax3 = subplot(3,1,3);
    S.ClickPlot = imagesc(S.ax3,squeeze(S.ProbMatrices(S.CodeNo,:,:))); colorbar;
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
        rectangle(S.ax3,'Position',S.gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
    end
    hold off
    
    S = getFigureTitle(S);    
    fprintf('yx=(%d,%d): code %d, %s\n', ...
        xy(2),xy(1),...
        S.CodeNo, S.GeneNames{S.CodeNo});
    
    if nCodesToUse>1
        S.sl = uicontrol('style','slide',...
            'unit','pix',...
            'position',[60 20 883 20],...
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

if strcmp(click_type,'normal')        
    LogProbPlot = log(conv(S.LambdaDist(:,S.CodeNo,b,r),S.HistProbs(:,b,r),'same'));
    PlotIdx = find(LogProbPlot>min(max(LogProbPlot)*5,-10));    %Only plot in range where prob above certain amount
    PlotIdx = min(PlotIdx):max(PlotIdx);    %So consecutive
    %Get background too
    HistZeroIndex = find(S.SymmHistValues == 0);
    BackgroundProb = log(S.HistProbs(HistZeroIndex+x,b,r));
    [~,MaxIdx] = max(BackgroundProb);
    BackgroundIdx1 = max(find(BackgroundProb==min(BackgroundProb)&find(BackgroundProb)<MaxIdx))+1;
    BackgroundIdx2 = min(find(BackgroundProb==min(BackgroundProb)&find(BackgroundProb)>MaxIdx))-1;
    BackgroundIdx = BackgroundIdx1:BackgroundIdx2;
    
    figure(35458);
    plot(x(BackgroundIdx),BackgroundProb(BackgroundIdx),'Color',[0.25, 0.25, 0.25],'LineWidth',0.1);
    
    %set(get(get(P1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on  
    plot(x(PlotIdx),LogProbPlot(PlotIdx),'Color',[0, 0.4470, 0.7410],'LineWidth',0.8);
    xline(f,'-','DisplayName','Spot Value','Color','red','LineWidth',1);   %Or xline(x(o.ZeroIndex-1+f))
    hold off
    leg1 = legend('$P(s\,\mid \,background)$','$P(s\,\mid \,gene\,\, and\,\, background)$','Spot Value');
    set(leg1,'Interpreter','latex');
    %legend('show');
    xlabel('Spot Intensity, $s$','interpreter','latex','FontSize',13)
    ylabel('Log Probability','interpreter','latex','FontSize',13);
    title('Probability distribution for '+string(S.GeneNames(S.CodeNo))+ ', round '+string(r)+' color channel '+string(b-1))
elseif strcmp(click_type,'alt')
    HistZeroIndex = find(S.SymmHistValues == 0);
    x2 = x(x<HistZeroIndex+f);      %So ensure indices>0
    hIndices = HistZeroIndex+f-x2;
    Use = hIndices<length(S.SymmHistValues);
    HistDist = S.HistProbs(hIndices(Use),b,r);
    LambdaIndices = find(x<HistZeroIndex+f);
    figure(9264992);
    plot(x(LambdaIndices(Use)),S.LambdaDist(LambdaIndices(Use),S.CodeNo,b,r));
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
S.ax2.Colorbar.Limits = [min(S.cBledCodes(S.CodeNo,:)),max(S.cBledCodes(S.CodeNo,:))];
S.ax2.Title.String = sprintf('Rank %d: Predicted Code for %s, code #%d',S.CodeIdx, S.GeneNames{S.CodeNo}, S.CodeNo);
S.ax2.YTick = 1:S.nBP;
S.ax2.YTickLabel = S.bpLabels;
S.ax2.YLabel.String = 'Color Channel';

    
S.ClickPlot = imagesc(S.ax3,squeeze(S.ProbMatrices(S.CodeNo,:,:)));
colorbar(S.ax3);
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

for ax = [S.ax1,S.ax2,S.ax3]
    hold on
    for r=1:S.nRounds
        rectangle(ax,'Position',S.gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':');
    end
    hold off
end
S = getFigureTitle(S); 
drawnow

end
