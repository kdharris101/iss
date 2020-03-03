function iss_change_plot(o,Method)
%Given issPlot3DObject, this function lets you change the details
%of the plot without closing the figure e.g. you can change
%o.CombiQualThresh or issPlot3DObject.ZThick to change the threshold value and the number
%of Z planes you see. 
%If Method == 'Prob', this changes the gene assignments to those given by
%the probability method rather than the dot product. In this case
%o.pScoreThresh is the threshold.

S = evalin('base', 'issPlot2DObject');
figure(S.FigNo);
h = findobj('type','line'); %KEY LINES: DELETE EXISTING SCATTER PLOTS SO CHANGE_SYMBOLS WORKS
delete(h);

if nargin<2 || isempty(Method)  
    if strcmpi(S.CallMethod,'DotProduct')
        S.QualOK = o.quality_threshold;
    elseif strcmpi(S.CallMethod,'Prob')
        S.QualOK = o.quality_threshold_prob;
    end
else
    if strcmpi('Prob',Method)
        S.SpotGeneName = o.GeneNames(o.pSpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        % which ones pass quality threshold (combi first)
        S.QualOK = o.quality_threshold_prob;
        S.CallMethod = 'Prob';
    else
        S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
        S.uGenes = unique(S.SpotGeneName);
        % which ones pass quality threshold (combi first)
        S.QualOK = o.quality_threshold;
        S.CallMethod = 'DotProduct';
    end
end
        
%S.SpotYXZ = o.SpotGlobalYXZ;
%S.Roi is the Roi for the current Z plane
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

assignin('base','issPlot2DObject',S)
