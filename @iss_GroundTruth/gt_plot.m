function gt_plot(o,Method)
%% o.gt_plot('Method');
% visually display ground truth spots found and missed.
%Sort gene numbers by round then channel.
S.GeneNoIndices = find(o.gtGeneNo>0);
[r,b] = find(o.gtGeneNo>0);
S.GeneNoIndices = sortrows([S.GeneNoIndices,r,b],[2,3]);
S.GeneNoIndices = S.GeneNoIndices(:,1);
gtImagesCell = cell(size(S.GeneNoIndices));
for i=1:length(S.GeneNoIndices)
    [r,b] = ind2sub(size(o.gtGeneNo),S.GeneNoIndices(i));
    gtImagesCell{i} = imread(o.gtBigImFiles{r,b})-o.TilePixelValueShift;
end
%Make dimension of all background images the same so zoom maintained.
XY_Dim = cell2mat(cellfun(@size,gtImagesCell,'uni',false));
gtImages = zeros([max(XY_Dim),length(S.GeneNoIndices)]);
for i=1:length(S.GeneNoIndices)
    gtImages(1:XY_Dim(i,1),1:1:XY_Dim(i,2),i) = gtImagesCell{i};
end

fh = figure(65828);
set(fh,'units','pixels','position',[500 200 800 600]);  %Left, Bottom, Width, Height
set(gcf, 'color', 'k');
set(gca, 'color', 'k');

S.Method = Method;
S.gtImages = gtImages;
S.i = 1;
S = rb_plot(o,S);

sl = uicontrol('style','slide',...
    'unit','pix',...
    'position',[60 8 693 18],...
    'min',1,'max',length(S.GeneNoIndices),'val',1,...
    'sliderstep',[1/(length(S.GeneNoIndices)-1) 1/(length(S.GeneNoIndices)-1)],...
    'callback',{@sl_call,o,S});
end

function S = rb_plot(o,S)
%Plots the round r, channel b ground truth as indicated by index S.i. 
ax = gca;
XLim = ax.XLim;
YLim = ax.YLim;
pf = o.CallMethodPrefix(S.Method);
S.Background = imagesc(S.gtImages(:,:,S.i));
colormap bone;
set(gca, 'YDir', 'normal');

[r,b] = ind2sub(size(o.gtGeneNo),S.GeneNoIndices(S.i));
caxis([0,3*o.gtColorTruePositiveThresh(r,b)]);
hold on
%Found ground truth
scatter(o.gtSpotGlobalYX{r,b}(o.([pf,'_gtFound']){r,b}==1,2),...
    o.gtSpotGlobalYX{r,b}(o.([pf,'_gtFound']){r,b}==1,1),30,'go','LineWidth',2);
%Missed ground truth
scatter(o.gtSpotGlobalYX{r,b}(o.([pf,'_gtFound']){r,b}==2,2),...
    o.gtSpotGlobalYX{r,b}(o.([pf,'_gtFound']){r,b}==2,1),30,'gx','LineWidth',2);
%pf method true positives
QualOK = quality_threshold(o,S.Method);
pfTP = o.([pf,'_gtIdentity']){r,b}==1 & o.([pf,'SpotCodeNo']) == o.gtGeneNo(r,b);
%Found
scatter(o.([pf,'SpotGlobalYX'])(pfTP&QualOK,2),o.([pf,'SpotGlobalYX'])(pfTP&QualOK,1),...
    20,[0 0.5 0],'o','LineWidth',1);
%Missed
scatter(o.([pf,'SpotGlobalYX'])(pfTP&~QualOK,2),o.([pf,'SpotGlobalYX'])(pfTP&~QualOK,1),...
    20,[0 0.5 0],'x','LineWidth',1);
%pf method false positives
pfFP = o.([pf,'_gtIdentity']){r,b}==2 & o.([pf,'SpotCodeNo']) == o.gtGeneNo(r,b) & QualOK;
scatter(o.([pf,'SpotGlobalYX'])(pfFP,2),o.([pf,'SpotGlobalYX'])(pfFP,1),...
    20,[0 0.5 0],'rx','LineWidth',1);
hold off
title(sprintf('Round %d, Channel %d: %s',r,b,o.GeneNames{o.gtGeneNo(r,b)}),'Color','white');
lgnd = legend('gtTP-Found','gtTP-Missed',[pf,'TP-PassQual'],[pf,'TP-FailQual'],[pf,'FP']);
set(lgnd,'color','black');
set(lgnd,'TextColor','white');
if XLim(2)~=1
    ax.XLim = XLim;
    ax.YLim = YLim;
end
end

function [] = sl_call(varargin)
% Callback for the slider.
[h,o,S] = varargin{[1,3,4]};
S.i = round(get(h,'value'));
S = rb_plot(o,S);
end
