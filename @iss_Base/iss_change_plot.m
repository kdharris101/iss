function iss_change_plot(o,Method,GeneType,GenesToShow,UseSpots)
%% iss_change_plot(o,Method,GeneType,GenesToShow,UseSpots)
%
% Given issPlot3DObject, this function lets you change the details
% of the plot without closing the figure e.g. you can change
% o.CombiQualThresh, o.pIntensityThresh, o.pScoreThresh, o.pScoreThresh2,
% o.pLogProbThresh
%
% o: iss object
% Method: 'DotProduct','Prob' or 'Pixel' to consider gene assignments given
% by o.SpotCodeNo, o.pSpotCodeNo and o.pxSpotCodeNo respectively.
% GeneType: Neuron or NonNeuron
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

if nargin<3 || isempty(GeneType)
    if ~isfield(S,'GeneType')
        GeneType = 'Neuron';
    else
        GeneType = S.GeneType;
    end
end
if ~strcmpi(GeneType,'Neuron') && ~strcmpi(GeneType,'NonNeuron')
    warning('Didn''t understand GeneType input so showing neuron type Genes');
    GeneType = 'Neuron';
end
S.GeneType = GeneType;


if nargin<4 || isempty(GenesToShow)
    GenesToShow = o.GeneNames;
    if ~isfield(S,'GeneNoToShow')
        %Only change if not previosuly given GenesToShow
        S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
    end
else
    S.GeneNoToShow = find(ismember(o.GeneNames,GenesToShow));
end


if nargin<2 || isempty(Method)  
    pf = o.CallMethodPrefix(S.CallMethod);
else
    if ~ismember({Method},o.CallMethods)
        error('Method invalid, must be member of o.CallMethods.');
    end
    S.CallMethod = Method;
    pf = o.CallMethodPrefix(S.CallMethod);
    S.SpotGeneName = o.GeneNames(o.([pf,'SpotCodeNo']));
    S.uGenes = unique(S.SpotGeneName);
end

if nargin>=5 && length(UseSpots)==length(o.([pf,'SpotCodeNo'])) && islogical(UseSpots)
    S.QualOK = UseSpots & ismember(o.([pf,'SpotCodeNo']),S.GeneNoToShow);
else
    if nargin>=5; warning('UseSpots not valid, using quality_threshold');end
    S.QualOK = quality_threshold(o,S.CallMethod) & ismember(o.([pf,'SpotCodeNo']),S.GeneNoToShow);
end

S.SpotYX = o.([pf,'SpotGlobalYX']);
        
InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));

%hold on; GET RID OF HOLD AND JUST DELETE PLOTS USING DELETE MEANS THAT THE
%ZOOM IS KEPT BETWEEN Z PLANES
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYX(MySpots,2), S.SpotYX(MySpots,1), '.', 'MarkerSize',...
            1,'Color',hsv2rgb([0 0 0.5]));
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    if strcmpi(GeneType,'Neuron')
        change_gene_symbols(0);
    elseif strcmpi(GeneType,'NonNeuron')
        change_gene_symbols_NonNeuron(0);
    end
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

assignin('base','issPlot2DObject',S)
