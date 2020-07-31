function iss_change_plot_contour(o,Method,GenesToShow,UseSpots)
%% iss_change_plot_contour(o,Method,GenesToShow,UseSpots)
%
% This function is identical to iss_change_plot(), except that it includes
% a semi-transparent contour plot overlayed onto the figure that shows the
% density of spots. Note that the color map is not absolute, the itensity
% of densities are not directly comparable between figures, although the
% relative spatial distribution should be.
%
% o: iss object
% Method: 'DotProduct','Prob' or 'Pixel' to consider gene assignments given
% by o.SpotCodeNo, o.pSpotCodeNo and o.pxSpotCodeNo respectively.
% GenesToShow: cell of gene names that you want to see e.g.
% [{'Npy'},{'Pvalb'}]. It is case sensitive.
% UseSpots: if you want to use your own thresholding, not
% o.quality_threshold. Logical array e.g. o.pSpotScore>0 or
% get_gene_clusters(o)

%%
S = evalin('base', 'issPlot2DObject');
figure(S.FigNo);
h = findobj('type','line'); %KEY LINES: DELETE EXISTING SCATTER PLOTS SO CHANGE_SYMBOLS WORKS
delete(h);


if nargin<3 || isempty(GenesToShow)
    GenesToShow = o.GeneNames;
    if ~isfield(S,'GeneNoToShow')
        %Only change if not previosuly given GenesToShow
        S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
    end
else
    S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
end


if nargin<2 || isempty(Method)  
    if strcmpi(S.CallMethod,'DotProduct')
        if nargin>=4 && length(UseSpots)==length(o.SpotScore) && islogical(UseSpots)
            S.QualOK = UseSpots & ismember(o.SpotCodeNo,S.GeneNoToShow);
        else
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.SpotCodeNo,S.GeneNoToShow);
        end
    elseif strcmpi(S.CallMethod,'Prob')
        if nargin>=4 && length(UseSpots)==length(o.pSpotScore) && islogical(UseSpots)
            S.QualOK = UseSpots & ismember(o.pSpotCodeNo,S.GeneNoToShow);    
        else
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pSpotCodeNo,S.GeneNoToShow);
        end
    elseif strcmpi(S.CallMethod,'Pixel')
        if nargin>=4 && length(UseSpots)==length(o.pxSpotScore) && islogical(UseSpots)
            S.QualOK = UseSpots & ismember(o.pxSpotCodeNo,S.GeneNoToShow);    
        else
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pxSpotCodeNo,S.GeneNoToShow);
        end
    end
else
    if strcmpi('Prob',Method)
        S.CallMethod = 'Prob';
        S.SpotGeneName = o.GeneNames(o.pSpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        if nargin>=4 && length(UseSpots)==length(o.pSpotScore) && max(UseSpots)==1
            S.QualOK = UseSpots & ismember(o.pSpotCodeNo,S.GeneNoToShow);  
        else
            % which ones pass quality threshold (combi first)
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold_prob');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pSpotCodeNo,S.GeneNoToShow);
        end
    elseif strcmpi('Pixel',Method)
        S.CallMethod = 'Pixel';
        S.SpotGeneName = o.GeneNames(o.pxSpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        if nargin>=4 && length(UseSpots)==length(o.pxSpotScore) && max(UseSpots)==1
            S.QualOK = UseSpots & ismember(o.pxSpotCodeNo,S.GeneNoToShow);  
        else
            % which ones pass quality threshold (combi first)
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold_prob');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.pxSpotCodeNo,S.GeneNoToShow);
        end
    else
        S.CallMethod = 'DotProduct';
        S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        if nargin>=4 && length(UseSpots)==length(o.SpotScore) && max(UseSpots)==1
            S.QualOK = UseSpots & ismember(o.SpotCodeNo,S.GeneNoToShow);    
        else
            % which ones pass quality threshold (combi first)
            if nargin>=4; warning('UseSpots not valid, using o.quality_threshold');end
            S.QualOK = o.quality_threshold(S.CallMethod) & ismember(o.SpotCodeNo,S.GeneNoToShow);
        end
    end
end

if strcmpi('Pixel',Method)
    S.SpotYX = o.pxSpotGlobalYX;
else
    S.SpotYX = o.SpotGlobalYX;
end
        
InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));

%hold on; GET RID OF HOLD AND JUST DELETE PLOTS USING DELETE MEANS THAT THE
%ZOOM IS KEPT BETWEEN Z PLANES
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYX(MySpots,2), S.SpotYX(MySpots,1), '.');
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

% all code for the contour plots are from here until the end
% the majority of which is adapted from content on Stack Overflow
% https://stackoverflow.com/questions/9134014/contour-plot-coloured-by-clustering-of-points-matlab
% Question asked by: HCAI https://stackoverflow.com/users/1134241/hcai
% Answer given by: Vidar https://stackoverflow.com/users/346645/vidar
SpotNo = find(ismember(o.GeneNames, GenesToShow));
Thresh = o.quality_threshold('Pixel') & ismember(o.pxSpotCodeNo, SpotNo);
Spots = o.pxSpotGlobalYX(Thresh,:);

Border = 0;
Sigma = 500;
StepSize = 100;

X = (Spots(:,2) / 1)';
Y = (Spots(:,1) / 1)';
D = [X' Y'];
N = length(X);

Xrange = [min(X)-Border max(X)+Border];
Yrange = [min(Y)-Border max(Y)+Border];

% set up coordinate grid
[XX, YY] = meshgrid(Xrange(1):StepSize:Xrange(2), Yrange(1):StepSize:Yrange(2));
YY = flipud(YY);

% parzen parameters and function handle
pf1 = @(C1,C2) (1/N)*(1/((2*pi)*Sigma^2)).* ...
         exp(-( (C1(1)-C2(1))^2+ (C1(2)-C2(2))^2)/(2*Sigma^2));

PPDF1 = zeros(size(XX));

% populate coordinate surface
[R, C] = size(PPDF1);
NN = length(D);
for c = 1:C
   for r = 1:R
       for d = 1:N
            PPDF1(r,c) = PPDF1(r,c) + ...
                pf1([XX(1,c) YY(r,1)], [D(d,1) D(d,2)]);
       end
   end
end

% normalize data
m1 = max(PPDF1(:));
PPDF1 = PPDF1 / m1;

% add PDF estimates to figure
OrigAx = gca;
ContourAx = axes;
[~, ContourH] = contourf(XX, YY, PPDF1);
colormap(gca, 'hot')
set(gca, 'YDir', 'reverse');
set(gca, 'XDir', 'reverse');
set(gca, 'color', 'none');
set(gca,'XTick',[], 'YTick', []);

drawnow; % Required for the below to work
FillsH = ContourH.FacePrims;
[FillsH.ColorType] = deal('truecoloralpha');
for i = 1:numel(FillsH)
   FillsH(i).ColorData(4) = 100;
end

set(gcf, 'currentaxes', OrigAx);

assignin('base','issPlot2DObject',S)
