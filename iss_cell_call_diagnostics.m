function iss_cell_call_diagnostics(o)
% iss_cell_call_diagnostics(o)
% run o.plot, then zoom into a region of interest, then run this.
%
% It will add cell class names for all cells. Then (once implemented) you
% can click on cells and it will give you diagnostics for them.

ProbThresh = 1e-3;
Tag = 'annotation';

% clear from last time
delete(findobj(gcf, 'tag', 'annotation'))


%% find spots in current axis frame
ax = axis; % gets current screen view in form [x0 x1 y0 y1]
xS = o.SpotGlobalYX(:,2); yS = o.SpotGlobalYX(:,1);
spotInFrame = (xS>=ax(1) & xS<=ax(2) & yS>=ax(3) & yS<= ax(4));
[topProb, SpotBestCell] = max(o.pSpotCell, [], 2);
% last "cell" is the noise cell, so exclude them as well as zero prob
spotOK = (SpotBestCell<size(o.pSpotCell,2) & topProb>0);
spotPlot = spotInFrame & spotOK;

spotCell = SpotBestCell(spotPlot);

h = plot([xS(spotPlot), o.CellYX(spotCell,2)]', ...
     [yS(spotPlot), o.CellYX(spotCell,1)]', ...
     'Color', [.3 .3 .3]);
 
set(h, 'tag', Tag);
        
uCells = unique(spotCell);
[~, BestClass] = max(o.pCellClass(uCells,:),[],2);            
h = text(o.CellYX(uCells,2), o.CellYX(uCells,1), o.ClassNames(BestClass),...
    'color', 'r', 'fontsize', 6);
set(h, 'tag', Tag);

%% now interactive mode

while(1)
    fprintf('Click on a cell, press key to quit\n')
    set(gca, 'Color', [1 1 1]*.2);
    [x,y,butt] = ginput(1);
    if butt~=1; break; end; % 
    set(gca, 'color', 'k');
    %    x = 5832; y = 7936;
    [~,myCell] = min(sum((o.CellYX-[y x]).^2,2));

    fprintf('------------------ Cell %d at %.0f,%.0f: -----------------\n', ...
            myCell, o.CellYX(myCell,1), o.CellYX(myCell,2));

    [uGenes,~,gInd] = unique(o.GeneNames);
    mySpots=find(o.pSpotCell(:,myCell)>ProbThresh);
    MyCellGeneCount = zeros(length(uGenes),1);
    for i=1:length(mySpots)
        ss = mySpots(i);
        MyProb = full(o.pSpotCell(ss,myCell));
        MyGeneNo = gInd(o.SpotCodeNo(ss));
        fprintf('Spot %d: %12s, with prob %.5f\n', ss, uGenes{MyGeneNo}, MyProb);
        MyCellGeneCount(MyGeneNo) = MyCellGeneCount(MyGeneNo)+ MyProb;
    end
    
    fprintf('-- Total Gene Count --\n');
    for i=1:length(MyCellGeneCount)
        if MyCellGeneCount(i)>1e-3
            fprintf('%12s:\t%f\n', uGenes{i}, MyCellGeneCount(i)); 
        end
    end
    fprintf('-- Class Posteriors --\n');
    for cc=find(o.pCellClass(myCell,:)>ProbThresh)
        fprintf('%20s:\t%.5f\n', o.ClassNames{cc}, o.pCellClass(myCell,cc)); 
    end
end
return
%%

MeanClassExp = zeros(nK, nG);
gSub = gSet.GeneSubset(GeneNames);
for k=1:nK-1 % don't include last since it is zero-expression class
    MeanClassExp(k,:) = o.Inefficiency * mean(gSub.ScaleCell(0).CellSubset(o.ClassNames{k}).GeneExp,2)';
end
lMeanClassExp = log(MeanClassExp + o.SpotReg); 

%%
nC = size(o.pSpotCell,2);
nN = o.nNeighbors+1; % last is misreads (always a neighbor)

    % CellGeneCount(nC, nG): number of copies of each gene in each cell
    CellGeneCount = zeros(nC,nG);
    for n=1:nN-1
        c = Neighbors(:,n);
        CellGeneCount = CellGeneCount + accumarray([c, SpotGeneNo], pSpotNeighb(:,n), [nC,nG]);
    end
    
    %% call cell gammas

    % eSpotGamma(nC, nK, nG); expected gamma parameter
    % elSpotGamma(nC, nK, nG); expected log gamma parameter
    ScaledMean = CellAreaFactor.*reshape(MeanClassExp,[1 nK nG]);
    eSpotGamma = (o.rSpot+reshape(CellGeneCount,[nC 1 nG]))./(o.rSpot + ScaledMean);
    elSpotGamma = psi(o.rSpot+reshape(CellGeneCount,[nC 1 nG])) - log(o.rSpot + ScaledMean); % expectation of log gamma


ScaledExp = reshape(MeanClassExp,[1 nK nG]) .* reshape(eGeneGamma,[1 1 nG]) .* CellAreaFactor   + o.SpotReg;

pNegBin = ScaledExp ./ (o.rSpot + ScaledExp);


figure(986543)
Myp = squeeze(pNegBin(myCell,:,:)); % nK by nG
WeightMap = MyCellGeneCount(:)' .* log(Myp) +  o.rSpot*log(1-Myp);
imagesc(WeightMap);
set(gca, 'xtick', 1:nG); set(gca, 'XTickLabel', GeneNames); set(gca, 'XTickLabelRotation', 90);
set(gca, 'ytick', 1:nK); set(gca, 'yTickLabel', o.ClassNames);
title(sprintf('Cell %d: Contribution of genes to class scores', myCell));

figure(19043765)
bar(eGeneGamma);
set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames);
set(gca, 'XTickLabelRotation', 90);
title('Gene efficiencies');
grid on

%%
DiagnosisPair = [];
if isempty(DiagnosisPair)
    [~, TopClasses] = sort(o.pCellClass(myCell,:), 'descend');
else
    TopClasses(1) = strmatch(DiagnosisPair{1}, o.ClassNames);
    TopClasses(2) = strmatch(DiagnosisPair{2}, o.ClassNames);;
end
while(1)
    GeneContrib = WeightMap(TopClasses(1),:) -  WeightMap(TopClasses(2),:);
    [sorted, order] = sort(GeneContrib);
    figure (986544);
    bar(sorted); 
    set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames(order));
    set(gca, 'XTickLabelRotation', 90);
    title(sprintf('Cell %d: Score for class %s vs %s', myCell, o.ClassNames{[TopClasses(1), TopClasses(2)]}));

    s = input('Enter a class to compare to, or return to quit: ', 's');
    if isempty(s); break; end
    TopClasses(2) = strmatch(s, o.ClassNames);

end