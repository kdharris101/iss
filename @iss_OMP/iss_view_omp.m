function SpotNo = iss_view_omp(o, FigNo, Norm, IncludeGT, SpotNum)
%% SpotNo = iss_view_omp2(o, FigNo, Norm, SpotNum)
%
% This function lets you view the spot code, predicted code and ompScore in
% each channel of UnbledCode of selected genes. Also shows you the omp
% coefficients for each gene and background code. 
%
% o: iss object
% FigNo: figure number (default, current figure)
% Norm: normalization option (described below)
% IncludeGT: if true, will also plot the ground truth rounds.
% SpotNum: number of spot to analyze or is yx location want to look
% at.(default is to click)
% SpotNo: returns the number of the spot analyzed.
%
% Norm = 1: Raw Colors
% Norm = 2: Z-score each round and colour channel. (SpotColors were in this
% form to find OMP coefficients).
% Norm = 3: Normalised by percentile for each color channel across all
% rounds


%%
pf = o.CallMethodPrefix('OMP');  

if nargin>=5
    if length(SpotNum)==1
        SpotNo = SpotNum;
    elseif length(SpotNum)==2
        xy = [SpotNum(2),SpotNum(1)];
        [~,SpotNo] = min(sum(abs(o.([pf,'SpotGlobalYX'])-SpotNum),2));
    end
    if nargin>=2 && ishandle(FigNo)
        figure(FigNo);
        xlim([o.([pf,'SpotGlobalYX'])(SpotNo,2)-20,o.([pf,'SpotGlobalYX'])(SpotNo,2)+20]);
        ylim([o.([pf,'SpotGlobalYX'])(SpotNo,1)-20,o.([pf,'SpotGlobalYX'])(SpotNo,1)+20]);
        rectangle('Position',[o.([pf,'SpotGlobalYX'])(SpotNo,2)-5,...
            o.([pf,'SpotGlobalYX'])(SpotNo,1)-5,10,10],'EdgeColor','g','LineWidth',4);
        drawnow;
    end
else
    if nargin>=2
        figure(FigNo);
    end
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    S = evalin('base', 'issPlot2DObject');
    InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
    PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
    [~,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);      
end

if nargin<4 || isempty(IncludeGT)
    IncludeGT = false;
end

%Different parameters for different methods
CodeNo = o.([pf,'SpotCodeNo'])(SpotNo);
%SpotColor is raw values (Norm 1)
SpotColor = o.([pf,'SpotColors'])(SpotNo,:,:);
CodeShape = size(SpotColor);
SpotGlobalYX = o.([pf,'SpotGlobalYX'])(SpotNo,:);
SpotCoefs = o.([pf,'Coefs'])(SpotNo,:);
SpotScore = o.([pf,'SpotScore'])(SpotNo);
SpotNeighbNonZero = o.([pf,'NeighbNonZeros'])(SpotNo);
SpotIntensity = o.([pf,'SpotIntensity'])(SpotNo);
%PredCode is scaled and shifted (Norm 2)
PredCode = SpotCoefs*o.([pf,'BledCodes'])(:,:);
PredCode = reshape(PredCode,CodeShape);

if nargin<3 || isempty(Norm)
    Norm = 2;
end

%Different Normalisations
if Norm == 1
    cSpotColor = double(SpotColor);
    cPredCode = PredCode.*o.z_scoreSCALE+o.z_scoreSHIFT;
elseif Norm == 2    
    cSpotColor = (double(SpotColor)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    cPredCode = PredCode;
elseif Norm == 3
    cSpotColor = double(SpotColor);
    cPredCode = PredCode.*o.z_scoreSCALE+o.z_scoreSHIFT;
    for b = 1:o.nBP
        bSpotColors = o.pSpotColors(:,b,:);
        p = double(prctile(bSpotColors(:), o.SpotNormPrctile));
        cSpotColor(:,b,:) = cSpotColor(:,b,:)/p;
        cPredCode(:,b,:) = cPredCode(:,b,:)/p;
    end
end

cSpotColor = squeeze(cSpotColor);
CodeShapeSqueeze = size(cSpotColor);
cPredCode = squeeze(cPredCode);
Error = cSpotColor - cPredCode;
ScoreMatrix = get_omp_score(o,(double(SpotColor)-o.z_scoreSHIFT)./o.z_scoreSCALE,...
    SpotCoefs,CodeNo);
v0 = min(-0.2,min(cSpotColor(:)));
v1 = max(0.2,max(cSpotColor(:)));
caxis_lims = [v0 v1];
score_caxis_lims = [-2, 2];

if IncludeGT
    gtSpotColor = o.([pf,'_gtColor'])(SpotNo,:,o.gtRounds);
    if Norm==3
        for b = 1:o.nBP
            bSpotColors = o.pSpotColors(:,b,:);
            p = double(prctile(bSpotColors(:), o.SpotNormPrctile));
            gtSpotColor(:,b,:) = gtSpotColor(:,b,:)/p;
        end
    elseif Norm==2
        %If normalising, z-score gt and then scale so has same max as
        %measured code.
        gtSHIFT = nanmedian(o.([pf,'_gtColor'])(:,:,o.gtRounds));
        gtSCALE = nanstd(o.([pf,'_gtColor'])(:,:,o.gtRounds));
        gtZSCORED = (o.([pf,'_gtColor'])(SpotNo,:,o.gtRounds)-gtSHIFT)./gtSCALE;
        gtSpotColor = squeeze(gtZSCORED);
    end
    MeasuredCode = [cSpotColor,squeeze(gtSpotColor)];
else
    MeasuredCode = cSpotColor;
end

%Get square outlining unbled code
gUnbled = reshape(o.UnbledCodes(CodeNo(1),:,:),CodeShapeSqueeze);
gSquares = zeros(o.nRounds,4);
for r=1:o.nRounds
    try
        gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
    end
end

try
    clf(430476599)
    figure(430476599)
catch
    figure(430476599)
end
subplot(4,1,1);
imagesc(MeasuredCode); colorbar('Color','w');
caxis(caxis_lims);
if Norm>1
    sq_color = 'g';
    colormap(gca,bluewhitered);
else
    sq_color = 'r';
end
title('Spot Code','Color','w');
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel','Color','w');
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor',sq_color,'LineWidth',2,'LineStyle',':')
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
                rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor',sq_color,'LineWidth',2,'LineStyle','-');
            else
                rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor','k','LineWidth',2,'LineStyle','-');
            end
        end
    end
end
hold off
set(gca,'XColor','w');
set(gca,'YColor','w');

subplot(4,1,2)
imagesc(cPredCode); colorbar('Color','w');
caxis(caxis_lims);
if Norm>1
    colormap(gca,bluewhitered);
end
title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo),'Color','w');
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel','Color','w');
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor',sq_color,'LineWidth',2,'LineStyle',':')
end
hold off
set(gca,'XColor','w');
set(gca,'YColor','w');

subplot(4,1,3);
imagesc(ScoreMatrix); colorbar('Color','w');
caxis(score_caxis_lims);
colormap(gca,bluewhitered);
title(sprintf('Score = %.2f', SpotScore),'Color','w');
set(gca, 'ytick', 1:o.nBP);
set(gca, 'YTickLabel', o.bpLabels);
ylabel('Color Channel','Color','w');
xlabel('Round','Color','w');
%title('$\ln\left({\frac{P(spot\,\mid \,gene\,\, and\,\, background)}{P(spot\,\mid \,background)}}\right)$','interpreter','latex','FontSize',15)
%title(sprintf('Log Probability the spot can be explained by gene - Log Probability it can be explained by background alone'));
%set(ClickPlot,'ButtonDownFcn',{@getCoord,o,SpotNo,CodeNo,SpotColor});
hold on
for r=1:o.nRounds
    rectangle('Position',gSquares(r,:),'EdgeColor','g','LineWidth',2,'LineStyle',':')
end
hold off
set(gca,'XColor','w');
set(gca,'YColor','w');

nCodes = size(o.CharCodes,1);
ClickPlot = subplot(4,1,4);
%scatter(1:nCodes,SpotCoefs(1:nCodes),46,[0, 0.4470, 0.7410],'.');
hold on
for i=1:nCodes
    S.h(i) = plot(i, SpotCoefs(i), '.');
end
legend(S.h,o.GeneNames);
legend off
change_gene_symbols(0);
%set(gcf, 'color', 'w');
set(gca, 'color', 'k');
h=findobj(gcf,'type','axes');
delete(h(1));
%scatter(InitialCodeNo,0.8,46,[0.8500, 0.3250, 0.0980],'.');
%scatter(CodeNo,SpotCoefs(CodeNo),66,[0.8500, 0.3250, 0.0980],'.');
scatter(nCodes+1:nCodes+o.nBackground,SpotCoefs(nCodes+1:nCodes+o.nBackground),46,'cyan','.');
%ylim([min([backWt-0.05,-0.05,gWt-0.05]),max([backWt+0.05,gWt+0.05,1])]);
ylabel('Coefficient','Color','w');
set(gca, 'xtick', [1:nCodes,nCodes+int8((o.nBackground+1)/2)]);
Labels = o.GeneNames;
Labels{end+1} = 'Background';
set(gca, 'XTickLabel', Labels);
set(gca,'XColor','w');
set(gca,'YColor','w');
ClickPlot.XTickLabel{CodeNo} = ['\color{red}',Labels{CodeNo}];
xtickangle(90);
hold off
set(ClickPlot,'ButtonDownFcn',{@getCoord,o,CodeShape,CodeShapeSqueeze,...
    sq_color,'g',ClickPlot.XTickLabel,SpotCoefs,cSpotColor,Norm,v0,v1,...
    cPredCode,SpotNo,CodeNo,score_caxis_lims,Error,pf,IncludeGT});

%Color different parameters depending if over threshold
% if SpotScore>o.pScoreThresh
%     c1 = [0,0.7,0]; else; c1 = [0,0,0];end
% if LogProbOverBackground<o.pLogProbThresh
%     c2 = [1,0,0]; else; c2 = [0,0,0];end
% if SpotScore+SpotScoreDev<o.pDevThresh
%     c3 = [1,0,0]; else; c3 = [0,0,0];end
% if SpotIntensity<o.pIntensityThresh
%     c4 = [1,0,0]; else; c4 = [0,0,0];end
% 
set(gcf,'Position',[409,92,793,707])
figtitle = sgtitle('', 'interpreter', 'tex','Color','w');   %'tex' required for colors
figtitle.String = sprintf('Spot %d: code %d, %s. Coef = %.2f, NeighbNonZero = %d, Score = %.2f, Intensity = %.0f',...
    SpotNo, CodeNo, o.GeneNames{CodeNo}, SpotCoefs(CodeNo), SpotNeighbNonZero,SpotScore,SpotIntensity);
%figtitle.Color='red';
%drawnow

fprintf('Spot %d at yx=(%d,%d): code %d, %s\n', ...
    SpotNo, SpotGlobalYX(1),SpotGlobalYX(2),...
    CodeNo, o.GeneNames{CodeNo});

if nargin>=4 && ishandle(FigNo)
    figure(FigNo);
    rect = findall(gcf,'Type', 'Rectangle');
    delete(rect);
end
end

function getCoord(aH,evnt,o,CodeShape,CodeShapeSqueeze,code_sq_color,error_sq_color,...
    original_xtick_labels,Coefs,cSpotColor,Norm,v0,v1,cPredCode,SpotNo,OrigCodeNo,...
    score_caxis_lims,Error,pf,IncludeGT)
%This plots the unbled codes
nCodes = size(o.CharCodes,1);
fig = ancestor(aH,'figure');
click_type = get(fig,'SelectionType');
ClickLoc = evnt.IntersectionPoint(1:2);
CodeNo = round(ClickLoc(1));
if CodeNo>nCodes
    CodeNo = nCodes;
elseif CodeNo < 1
    CodeNo = 1;
end

rect = findall(gcf,'Type', 'Rectangle');
caxis_lims = [v0 v1];
error_caxis_lims = [min(v0,min(Error(:))), max(v1,max(Error(:)))];
h=findobj(gcf,'type','axes'); 

if strcmp(click_type,'normal')
        
    set(gcf,'CurrentAxes',h(5-2));
    imagesc(cPredCode); colorbar('Color','w'); 
    caxis(caxis_lims);
    if Norm>1
        colormap(gca,bluewhitered);
    end
    title(sprintf('Predicted Code for %s, code #%d. Coef = %.2f', o.GeneNames{CodeNo}, CodeNo, Coefs(CodeNo)),'Color','w');     
    
    if CodeNo ~= OrigCodeNo
        if Coefs(CodeNo)==0
            ScoreMatrix = get_omp_score(o,(double(o.([pf,'SpotColors'])(SpotNo,:,:))-o.z_scoreSHIFT)./o.z_scoreSCALE,...
                o.([pf,'Coefs'])(SpotNo,:),OrigCodeNo);
            CodeNoPrint = OrigCodeNo;
        else
            ScoreMatrix = get_omp_score(o,(double(o.([pf,'SpotColors'])(SpotNo,:,:))-o.z_scoreSHIFT)./o.z_scoreSCALE,...
                o.([pf,'Coefs'])(SpotNo,:),CodeNo);
            CodeNoPrint = CodeNo;
        end
        SpotScore = sum(ScoreMatrix(:));
        set(gcf,'CurrentAxes',h(5-3));
        imagesc(ScoreMatrix); colorbar('Color','w');
        caxis(score_caxis_lims);
        colormap(gca,bluewhitered);     
        title(sprintf('Gene %d, %s Score = %.2f', CodeNoPrint, o.GeneNames{CodeNoPrint}, SpotScore),'Color','w');
    else

        set(gcf,'CurrentAxes',h(5-3));
        imagesc(Error); colorbar('Color','w');
        caxis(error_caxis_lims);
        colormap(gca,bluewhitered);
        title(sprintf('Error with %s, code #%d, coef = %.2f is %.2f', o.GeneNames{CodeNo}, CodeNo, Coefs(CodeNo),sum(abs(Error(:)))),'Color','w');        
    end
   
    
    
elseif strcmp(click_type,'alt') 
    
    %Change PredCode and Error to without gene clicked on
    SpotCoefsClick = Coefs;
    SpotCoefsClick(CodeNo) = 0;
    PredCodeClick = SpotCoefsClick*o.([pf,'BledCodes'])(:,:);
    PredCodeClick = reshape(PredCodeClick,CodeShape);
    %Different Normalisations
    if Norm == 1
        cPredCodeClick = PredCodeClick.*o.z_scoreSCALE+o.z_scoreSHIFT;
    elseif Norm == 2
        cPredCodeClick = PredCodeClick;
    elseif Norm == 3
        cPredCodeClick = PredCodeClick.*o.z_scoreSCALE+o.z_scoreSHIFT;
        for b = 1:o.nBP
            bSpotColors = o.([pf,'SpotColors'])(:,b,:);
            p = double(prctile(bSpotColors(:), o.SpotNormPrctile));
            cPredCodeClick(:,b,:) = cPredCodeClick(:,b,:)/p;
        end
    end
    cPredCodeClick = squeeze(cPredCodeClick);
    ErrorClick = cSpotColor - squeeze(cPredCodeClick);
    
    set(gcf,'CurrentAxes',h(5-2));
    imagesc(cPredCodeClick); colorbar('Color','w'); 
    caxis(caxis_lims);
    if Norm>1
        colormap(gca,bluewhitered);
    end
    title(sprintf('Predicted Code without %s, code #%d. Coef = %.2f',...
        o.GeneNames{CodeNo}, CodeNo, Coefs(CodeNo)),'Color','w');  
    
    set(gcf,'CurrentAxes',h(5-3));
    imagesc(ErrorClick); colorbar('Color','w'); 
    caxis(error_caxis_lims);
    colormap(gca,bluewhitered);
    title(sprintf('Error without %s, code #%d, coef = %.2f is %.2f',...
        o.GeneNames{CodeNo}, CodeNo, Coefs(CodeNo),sum(abs(ErrorClick(:)))),'Color','w');
      
             
end
    %Delete current rectangles    
    delete(rect);
    
    %Get square outlining unbled code
    gUnbled = reshape(o.UnbledCodes(CodeNo,:),CodeShapeSqueeze);
    gSquares = zeros(o.nRounds,4);
    for r=1:o.nRounds
        try
            gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
        end
    end
        
    for subplot = 1:3
        set(gcf,'CurrentAxes',h(5-subplot));
        if subplot==3
            sq_color = error_sq_color; else; sq_color = code_sq_color;
        end
        hold on
        for r=1:o.nRounds
            rectangle('Position',gSquares(r,:),'EdgeColor',sq_color,'LineWidth',2,'LineStyle',':')
        end
        if IncludeGT && subplot==1
            for r=o.gtRounds
                for b=1:o.nBP
                    if o.gtGeneNo(r,b)==0; continue; end
                    if o.gtGeneNo(r,b)==CodeNo
                        rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor',sq_color,'LineWidth',2,'LineStyle','-');
                    else
                        rectangle('Position',[r-0.5,b-0.5,1,1],'EdgeColor','k','LineWidth',2,'LineStyle','-');
                    end
                end
            end
        end
        hold off
        
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'YTickLabel', o.bpLabels);
        ylabel('Color Channel','Color','w');
        set(gca,'XColor','w');
        set(gca,'YColor','w');         
    end  
    h(1).XTickLabel = original_xtick_labels;
    h(1).XTickLabel{CodeNo} = ['\color{green}',o.GeneNames{CodeNo}];
    drawnow;  
end

