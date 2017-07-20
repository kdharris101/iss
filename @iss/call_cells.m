function o = call_cells(o, gSet)
% o = o.call_cells(gSet)
%  
% Cell calling via negative binomial model
% 
% input gSet: GeneSet structure containing results of scRNAseq clustering
%
% Creates outputs: 
% pCellClass: posterior probability of each cell to be in each class (nCells*nClasses)
% pSpotCell: posterior probability of each spot to be in top 5 neighboring cells (nSpots * nCells, sparse)
% note that the last class is a zero-expressing cell; the last cell is background

%% diagnostic parameters
o.CellCallShowCenter = [1620 1650];
o.ExampleCellCenter = o.CellCallShowCenter;
Class1 = 'Calb2.Vip.Igfbp4';
Class2 = 'Zero';

%% exclude genes that are useless - and that's none of them?
ExcludeGenes = {}; %{'Vsnl1', 'Atp1b1', 'Slc24a2', 'Tmsb10'};

%% load properties of region of interest
load(o.CellMapFile); % CellMap and y0, y1, x0, x1 that are its coords in full image

IncludeSpot = ~ismember(SpotGenes, ExcludeGenes) ...
    & inpolygon(o.SpotGlobalYX(:,1), o.SpotGlobalYX(:,2), o.CellCallRegionYX(:,1), o.CellCallRegionYX(:,2));
SpotYX = SpotYX(IncludeSpot,:) - [y0 x0] + 1;
SpotGenes = SpotGenes(IncludeSpot);


%% get info about cells
rp = regionprops(CellMap);
CellYX = fliplr(vertcat(rp.Centroid)); % convert XY to YX
CellArea0 = vertcat(rp.Area); 

MeanCellRadius = mean(sqrt(CellArea0/pi))*.5; % the dapi part is only half of the typical radius
RelCellRadius = [sqrt(CellArea0/pi)/MeanCellRadius; 1]; % but here we want the whole thing
%RelCellArea = [ones(size(CellArea0)) ; 1]; % last "cell" is background signal 

%% get arrays ready

% SpotGene(nS): which gene is each spot
% MeanClassExp(nK,nG): mean expression of each gene in each class
% Neighbors(nS, nN): closest neighboring cells for each spot
% D(nS, nN): distance penalty for each of these
% GeneNames(nG): name of each gene
% ClassNames(nK): name of each class

[GeneNames, ~, SpotGene] = unique(SpotGenes);
TotGeneSpots = accumarray(SpotGene,1);
[ClassNames] = vertcat(unique(gSet.Class, 'stable'), {'Zero'});

nG = length(GeneNames);
nK = length(ClassNames); % last is zero-expression
nC = size(CellYX,1)+1; % last is misreads
nS = size(SpotYX,1);
nN = o.nNeighbors+1; % last is misreads (always a neighbor)

o.ClassPrior = [.5*ones(1,nK-1)/nK .5];

ClassDisplayNames = ClassNames;

MeanClassExp = zeros(nK, nG);
gSub = gSet.GeneSubset(GeneNames);
for k=1:nK-1 % don't include last since it is zero-expression class
    MeanClassExp(k,:) = o.Inefficiency * mean(gSub.ScaleCell(0).CellSubset(ClassNames{k}).GeneExp,2)';
end
% MeanClassExp = MeanClassExp;
lMeanClassExp = log(MeanClassExp + o.SpotReg); 

% now find each spot's neighboring cells and distances
[Neighbors, Dist] = knnsearch(CellYX, SpotYX, 'K', nN);
Neighbors(:,end) = nC; % set last neighbor to misreads

D = -Dist.^2./(2*MeanCellRadius^2) - log(2*pi*MeanCellRadius^2); % don't normalize: bigger cells express more
D(:,end) = log(o.MisreadDensity); % this is log likelihood of misread

% any inside cell radius given a bonus
SpotInCell = IndexArrayNan(CellMap, SpotYX');
if Neighbors(SpotInCell>0,1)~=SpotInCell(SpotInCell>0)
    error('a spot is in a cell not closest neighbor!');
end
D(SpotInCell>0, 1) = D(SpotInCell>0, 1) + o.InsideCellBonus;

LogClassPrior = log(o.ClassPrior);

% this is area factor relative to that of the average cell
CellAreaFactor = (exp(-RelCellRadius.^2/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus)) ...
    / (exp(-1/2)*(1-exp(o.InsideCellBonus)) + exp(o.InsideCellBonus));


%% initialize variables for main loop

pSpotNeighb = zeros(nS, nN); % prob each spot goes to each neighboring cell: last assigned to noise
pCellClass = zeros(nC, nK); % prob each cell goes to each class: last has zero expression
% initialize arrays

% start a spot in cell it is in, otherwise misread
pSpotNeighb(Neighbors==SpotInCell)=1;
pSpotNeighb(SpotInCell==0,end)=1;

% gammas start off as priors
eSpotGamma = ones(nC, nK, nG);
elSpotGamma = ones(nC, nK, nG)*psi(1); % start with r=1 prior, no spots

eGeneGamma = ones(nG,1); % start with just 1

% this is to check convergence
pSpotNeighbOld = zeros(nS, nN);

for i=1:o.CellCallMaxIter
    % CellGeneCount(nC, nG): number of copies of each gene in each cell
    CellGeneCount = zeros(nC,nG);
    for n=1:nN-1
        c = Neighbors(:,n);
        CellGeneCount = CellGeneCount + accumarray([c, SpotGene], pSpotNeighb(:,n), [nC,nG]);
    end

    %% call cells

    
    % ScaledExp(nC, nK, nG): expected expression under current parameters
    ScaledExp = reshape(MeanClassExp,[1 nK nG]) .* reshape(eGeneGamma,[1 1 nG]) .* CellAreaFactor   + o.SpotReg;

    % pNegBin(nC, nK, nG): negbin parameter
    pNegBin = ScaledExp ./ (o.rSpot + ScaledExp);
    
    % wCellClass(nC, nK): summed log likelihoods
    wCellClass = sum(reshape(CellGeneCount,[nC 1 nG]).*log(pNegBin) + o.rSpot*log(1-pNegBin),3) + LogClassPrior;
    
    % pCellClass(nC, nK): probabilities
    pCellClass = LogLtoP(wCellClass')';

    
    %% call cell gammas

    % eSpotGamma(nC, nK, nG); expected gamma parameter
    % elSpotGamma(nC, nK, nG); expected log gamma parameter
    ScaledMean = CellAreaFactor.*reshape(MeanClassExp,[1 nK nG]);
    eSpotGamma = (o.rSpot+reshape(CellGeneCount,[nC 1 nG]))./(o.rSpot + ScaledMean);
    elSpotGamma = psi(o.rSpot+reshape(CellGeneCount,[nC 1 nG])) - log(o.rSpot + ScaledMean); % expectation of log gamma

    %% call gene gammas (etas)
    % to count non-background expression of each gene
    BackgroundSpots = accumarray(SpotGene, pSpotNeighb(:,end), [nG 1]);
    % total predicted by other models
    TotPredicted = sum(shiftdim(sum(eSpotGamma.*pCellClass.*CellAreaFactor,1),1).*MeanClassExp,1)';
    eGeneGamma = (o.rGene + TotGeneSpots - BackgroundSpots)./(o.rGene + TotPredicted);
    if 0
        for gg=1:nG; 
            fprintf('%s:\t%f\n', GeneNames{gg}, eGeneGamma(gg)); 
        end
    end
    
    %% call spots
    % wSpotCell(nS, nN)
    aSpotCell = zeros(nS, nN);
    for n=1:nN-1 % don't include misread possibility
        c = Neighbors(:,n);
        aSpotCell(:,n) = sum(pCellClass(c,:) .* lMeanClassExp(:,SpotGene)',2) + ...
            sum(pCellClass(c,:) .* bi(elSpotGamma, c, 1:nK, SpotGene), 2);
    end
    wSpotCell = aSpotCell + D ;
    
    pSpotNeighb = LogLtoP(wSpotCell')';
    MeanProbChanged = max(abs(pSpotNeighb(:)-pSpotNeighbOld(:)));
    fprintf('Iteration %d, mean prob change %f\n', i, MeanProbChanged)
    Converged = ( MeanProbChanged<o.CellCallTolerance);
    pSpotNeighbOld = pSpotNeighb;
    
    
    %% diagnostics
    if ~isempty(o.CellCallShowCenter) && (Converged || o.Graphics==2);
        figure(3985471); cla
        iss_make_figure(o, SpotYX, SpotGenes); 
        hold on
        [~, BestNeighb] = max(pSpotNeighb,[],2);
        %BestNeighb = ones(nS,1);
        SpotBestNeighb = bi(Neighbors,(1:nS)',BestNeighb(:));
        rn = SpotBestNeighb<nC & sum(abs(SpotYX-o.CellCallShowCenter).^2,2)<o.CellCallShowRad.^2;
        plot([SpotYX(rn,2) , CellYX(SpotBestNeighb(rn),2)]', ...
            [SpotYX(rn,1) , CellYX(SpotBestNeighb(rn),1)]', 'Color', [.3 .3 .3]);
        [~, BestClass] = max(pCellClass(1:end-1,:),[],2);            
        text(CellYX(:,2), CellYX(:,1), ClassDisplayNames(BestClass), 'color', 'r', 'fontsize', 6);
    %     input('press key', 's');
        axis(reshape(fliplr([o.CellCallShowCenter;o.CellCallShowCenter])+[-o.CellCallShowRad; o.CellCallShowRad], 1, 4));

        if ~isempty(o.ExampleCellCenter)
            [~, MyCell] = min(sum((CellYX-o.ExampleCellCenter).^2,2));
            fprintf('------------------ Cell %d: -----------------\n', MyCell);
            for ss=find((Neighbors(:,1)==MyCell))'
                fprintf('Spot %d: %s, with prob %f\n', ss, GeneNames{SpotGene(ss)}, pSpotNeighb(ss,1));
            end
            fprintf('-- Total Gene Count --\n');
            for gg=find(CellGeneCount(MyCell,:)>1e-3) 
                fprintf('%s:\t%f\n', GeneNames{gg}, CellGeneCount(MyCell,gg)); 
            end
            fprintf('-- Class Posteriors --\n');
            for cc=find(pCellClass(MyCell,:)>1e-3)
                fprintf('%s:\t%e\n', ClassDisplayNames{cc}, pCellClass(MyCell,cc)); 
            end
            
            figure(986543)
            Myp = squeeze(pNegBin(MyCell,:,:)); % nK by nG
            WeightMap = CellGeneCount(MyCell,:) .* log(Myp) +  o.rSpot*log(1-Myp);
            imagesc(WeightMap);
            set(gca, 'xtick', 1:nG); set(gca, 'XTickLabel', GeneNames); set(gca, 'XTickLabelRotation', 90);
            set(gca, 'ytick', 1:nK); set(gca, 'yTickLabel', ClassNames);
            title(sprintf('Cell %d: Contribution of genes to class scores', MyCell));
            
            figure(19043765)
            bar(eGeneGamma);
            set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames);
            set(gca, 'XTickLabelRotation', 90);
            title('Gene efficiencies');
            grid on

            
            if 1
                GeneContrib = WeightMap(strmatch(Class1, ClassNames),:) ...
                    - WeightMap(strmatch(Class2, ClassNames),:);
                [sorted order] = sort(GeneContrib);
                figure (986544);
                bar(sorted); 
                set(gca, 'XTick', 1:nG), set(gca, 'XTickLabel', GeneNames(order));
                set(gca, 'XTickLabelRotation', 90);
                title(sprintf('Cell %d: Score for class %s vs %s', MyCell, Class1, Class2));
            end

        end
%         keyboard
    end
    
    %%
    if Converged; break; end
end

%% make dense array output

o.pSpotCell = sparse(repmat(1:nS,1,nN)', Neighbors(:), pSpotNeighb(:));
