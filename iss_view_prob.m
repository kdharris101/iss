function SpotNo = iss_view_prob(o, FigNo, Norm, Method, SpotNum)
    %This function lets you view the spot code, gene code,
    %spotcode-lambda*bledcode and ln(Prob) for the best match for a chosen
    %spot.
    %Now have option for different normalisations
    %Norm = 1: Raw Colors
    %Norm = 2: Normalised by o.SpotNormPrctile in each colour channel and round,
    %then if o.CallSpotsCodeNorm='WholeCode', normalise so whole code has L2 norm
    %of 1 but if o.CallSpotsCodeNorm='Round', normalise so each round has L2 norm of 1.
    %Norm = 3: Normalised by percentile for each color channel across all
    %rounds
    %Method = 'Prob' or 'Pixel' to consider gene assignments given
    %by o.pSpotCodeNo or o.pxSpotCodeNo respectively.
    
    
    if nargin>=5
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
            if strcmpi('Prob',Method)
                S.SpotYX = o.SpotGlobalYX;
            elseif strcmpi('Pixel',Method)
                S.SpotYX = o.pxSpotGlobalYX;
            end
            S.QualOK = 1;
        end
        InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
        PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
        [~,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
        SpotNo = PlotSpots(SpotIdx);        
            

    end
    
    %Different parameters for different methods
    if strcmpi('Prob',Method)
        CodeNo = o.pSpotCodeNo(SpotNo);
        SpotColor = o.cSpotColors(SpotNo,:,:);
        SpotScore = o.pSpotScore(SpotNo);
        LogProbOverBackground = o.pLogProbOverBackground(SpotNo);
        SpotScoreDev = o.pSpotScoreDev(SpotNo);
        SpotIntensity = o.pSpotIntensity(SpotNo); 
        SpotGlobalYX = o.SpotGlobalYX(SpotNo,:);
    elseif strcmpi('Pixel',Method)
        CodeNo = o.pxSpotCodeNo(SpotNo);
        SpotColor = o.pxSpotColors(SpotNo,:,:);
        SpotScore = o.pxSpotScore(SpotNo);
        LogProbOverBackground = o.pxLogProbOverBackground(SpotNo);
        SpotScoreDev = o.pxSpotScoreDev(SpotNo);
        SpotIntensity = o.pxSpotIntensity(SpotNo);
        SpotGlobalYX = o.pxSpotGlobalYX(SpotNo,:);
    else
        error('Spot calling method not valid, should be Prob or Pixel');
    end
    
    if nargin<3 || isempty(Norm)
        Norm = 1;
    end
    
    %Different Normalisations
    if isempty(Norm) || Norm == 1
        cSpotColor = SpotColor;
        cBledCodes = o.pBledCodes;
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
        cBledCodes = bsxfun(@rdivide, o.BledCodes, sqrt(sum(o.BledCodes.^2,2)));
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
        cBledCodes = change_bled_codes(o,NewBleedMatrix);
    end

    MeasuredCode = squeeze(cSpotColor);
    CodeShape = size(MeasuredCode);
    BledCode = cBledCodes(CodeNo,:);
    ProbMatrix = get_prob_matrix(o,squeeze(SpotColor),CodeNo);
    
    try
        clf(430476533)
        figure(430476533)
    catch
        figure(430476533)
    end
    subplot(3,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Spot Code'));
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    
    subplot(3,1,2)
    imagesc(reshape(BledCode, CodeShape)); colorbar
    %caxis([0 max(cBledCode(:))]);
    title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    
    ClickPlot = subplot(3,1,3);
    ClickPlot(1) = imagesc(ProbMatrix); colorbar
    %caxis([min(ProbMatrix(:)) max(ProbMatrix(:))]);
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');
   	title('$\ln\left({\frac{P(spot\,\mid \,gene\,\, and\,\, background)}{P(spot\,\mid \,background)}}\right)$','interpreter','latex','FontSize',15)
    %title(sprintf('Log Probability the spot can be explained by gene - Log Probability it can be explained by background alone'));
    set(ClickPlot,'ButtonDownFcn',{@getCoord,o,SpotNo,CodeNo,SpotColor});
    
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
x = min(f,min(o.cSpotColors(:))-1):max(f,max(o.cSpotColors(:))-1);

if strcmp(click_type,'normal')        
    LogProbPlot = log(conv(o.LambdaDist(:,CodeNo,b,r),o.HistProbs(:,b,r),'same'));
    PlotIdx = find(LogProbPlot>min(max(LogProbPlot)*5,-10));    %Only plot in range where prob above certain amount
    PlotIdx = min(PlotIdx):max(PlotIdx);    %So consecutive
    %Get background too
    HistZeroIndex = find(o.SymmHistValues == 0);
    BackgroundProb = log(o.HistProbs(HistZeroIndex+x,b,r));
    [~,MaxIdx] = max(BackgroundProb);
    BackgroundIdx1 = max(find(BackgroundProb==min(BackgroundProb)&find(BackgroundProb)<MaxIdx))+1;
    BackgroundIdx2 = min(find(BackgroundProb==min(BackgroundProb)&find(BackgroundProb)>MaxIdx))-1;
    BackgroundIdx = BackgroundIdx1:BackgroundIdx2;
    
    figure(35428);
    plot(x(BackgroundIdx),BackgroundProb(BackgroundIdx),'Color',[0.25, 0.25, 0.25],'LineWidth',0.1);
    
    %set(get(get(P1(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on  
    plot(x(PlotIdx),LogProbPlot(PlotIdx),'Color',[0, 0.4470, 0.7410],'LineWidth',0.8);
    xline(f,'-','DisplayName','Spot '+string(SpotNo)+ ' Value','Color','red','LineWidth',1);   %Or xline(x(o.ZeroIndex-1+f))
    hold off
    leg1 = legend('$P(s\,\mid \,background)$','$P(s\,\mid \,gene\,\, and\,\, background)$','Spot '+string(SpotNo)+ ' Value');
    set(leg1,'Interpreter','latex');
    %legend('show');
    xlabel('Spot Intensity, $s$','interpreter','latex','FontSize',13)
    ylabel('Log Probability','interpreter','latex','FontSize',13);
    title('Probability distribution for '+string(o.GeneNames(CodeNo))+ ', round '+string(r)+' color channel '+string(b-1))
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
    
end
end
