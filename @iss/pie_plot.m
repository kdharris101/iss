function iss_pie_plot(o, CellYX, pCellClass, ClassNames, Colors)
% plot pie chart for each cell, showing probability of it belonging to all
% classes

o.MinPieProb = .1;
o.PieSize = 5;

nC = size(CellYX,1);
nK = size(pCellClass,2);

% Colors0 = hsv(ceil(nK*1.2));
% Colors = [Colors0(1:nK-1,:); 0 0 0];

clf; 
set(gcf, 'Color', 'k');
set(gca, 'color', 'k');
hold on
for c=1:nC
    

    pMy = pCellClass(c,:);
    WorthShowing = find(pMy>o.MinPieProb);

    h = pie(pMy(WorthShowing), repmat({''}, sum(WorthShowing>0)));

    for i=1:length(h)/2
        hno = (i*2-1);
%         Index = find(strcmp(h(i*2).String, NickNames), 1);
        set(h(hno), 'FaceColor', Colors(WorthShowing(i),:));
        set(h(hno), 'Xdata', get(h(hno), 'Xdata')*o.PieSize + CellYX(c,2));
        set(h(hno), 'Ydata', get(h(hno), 'Ydata')*o.PieSize + CellYX(c,1));
        set(h(hno), 'EdgeAlpha', 0);
    end
    
    if mod(c,2000)==0
        drawnow
    end
end

yMax = max(CellYX(:,1));
xMax = max(CellYX(:,2));
yMin = min(CellYX(:,1));
xMin = min(CellYX(:,2));

if nargin>=4
    ClassShown = find(any(pCellClass>o.MinPieProb,1));
    nShown = length(ClassShown);
    for k=1:nShown
        h = text(xMax*1.1 - xMin*.1, yMin + k*(yMax-yMin)/nShown, ClassNames{ClassShown(k)});
        set(h, 'color', Colors(ClassShown(k),:));
    end
end

