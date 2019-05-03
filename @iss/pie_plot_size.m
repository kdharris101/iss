function pie_plot(o)
% plot pie chart for each cell, showing probability of it belonging to all
% classes

nC = size(o.CellYX,1);
nK = size(o.pCellClass,2);

% find classes to collapse
CollapseMe = zeros(nK,1);
Colors = zeros(nK,3);
DisplayName = o.ClassNames;
for i=1:size(o.ClassCollapse,1)
    ClassList = o.ClassCollapse{i,1};
    for j=1:length(ClassList)
        MyClasses = strmatch(ClassList{j}, o.ClassNames);
        if length(MyClasses)==0; continue; end
        CollapseMe(MyClasses)=i;
        Colors(MyClasses,:) = repmat(o.ClassCollapse{i,3},length(MyClasses),1);
        DisplayName(MyClasses) = o.ClassCollapse(i,2);
    end
end

nColorWheel = sum(CollapseMe==0);

Colors0 = hsv(ceil(nColorWheel*1.2));
Colors(~CollapseMe,:) = Colors0(1:nColorWheel,:); % last is zero

figure(43908765)
% figure
clf; 
set(gcf, 'Color', 'w');
set(gca, 'color', 'w');
hold on

% load(o.CellMapFile, 'RelCellRadius');

for c=1:nC
   

    pMy = o.pCellClass(c,:);
    WorthShowing = find(pMy>o.MinPieProb);
    if ~isempty(WorthShowing)

        h = pie(pMy(WorthShowing), repmat({''}, sum(WorthShowing>0)));

        for i=1:length(h)/2
            hno = (i*2-1);
    %         Index = find(strcmp(h(i*2).String, NickNames), 1);
            set(h(hno), 'FaceColor', Colors(WorthShowing(i),:));
            
%             % all same size
%             set(h(hno), 'Xdata', get(h(hno), 'Xdata')*o.PieSize + o.CellYX(c,2));
%             set(h(hno), 'Ydata', get(h(hno), 'Ydata')*o.PieSize + o.CellYX(c,1));

%             % size based on relative cell radius
%             set(h(hno), 'Xdata', get(h(hno), 'Xdata')*o.PieSize*RelCellRadius(c) + o.CellYX(c,2));
%             set(h(hno), 'Ydata', get(h(hno), 'Ydata')*o.PieSize*RelCellRadius(c) + o.CellYX(c,1));            

            % size based on number of reads
            set(h(hno), 'Xdata', get(h(hno), 'Xdata')*o.PieSize*sum(o.pSpotCell(:,c)) + o.CellYX(c,2));
            set(h(hno), 'Ydata', get(h(hno), 'Ydata')*o.PieSize*sum(o.pSpotCell(:,c)) + o.CellYX(c,1));            
            
%             set(h(hno), 'EdgeAlpha', 0);
            set(h(hno), 'EdgeAlpha', 1, 'LineWidth', .1, 'EdgeColor', [.6 .6 .6]);
        end
    end
    
    if mod(c,2000)==0
        drawnow
    end
end

yMax = max(o.CellYX(:,1));
xMax = max(o.CellYX(:,2));
yMin = min(o.CellYX(:,1));
xMin = min(o.CellYX(:,2));

ClassShown = find(any(o.pCellClass>o.MinPieProb,1));
ClassDisplayNameShown = DisplayName(ClassShown);
[uDisplayNames, idx] = unique(ClassDisplayNameShown, 'stable');
nShown = length(uDisplayNames);
for k=1:nShown
    h = text(xMax*1.1 - xMin*.11, yMin + k*(yMax-yMin)/nShown, DisplayName{ClassShown(idx(k))}, 'fontsize', 8);
    set(h, 'color', Colors(ClassShown(idx(k)),:));
end

