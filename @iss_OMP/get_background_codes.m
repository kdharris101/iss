function [BackgroundEigenvectors,BackgroundEigenvalues,BackgroundMaxGeneDotProduct,...
    BackgroundMaxGeneDotProductGene,EigenvectorTiles] = get_background_codes(o,ChannelStrips,DotProductThresh,nTilesToUse)
%% Find vectors to best represent background
%First find all pixels not well represented by genes. Then find
%eigenvectors of covariance matrix of resulting pixels - these are the
%background eigenvectors.
%o: iss_OMP object.
%ChannelStrips: if true, returns background eigenvectors that are just a
%strip in each channel. 
%DotProductThresh: If best dot product less than this, then is background
%nTilesToUse: number of tiles to consider to find background eigenvectors. 
%EigenvectorTiles: tiles used to find the background eigenvectors
%%
if nargin<2 || isempty(ChannelStrips)
    ChannelStrips = false;
end
nCodes = length(o.CharCodes);
NormBledCodes = o.z_scoreBledCodes./vecnorm(o.z_scoreBledCodes(:,:),2,2);

if ChannelStrips
    BackgroundEigenvectors = zeros(o.nBP,o.nBP,o.nRounds);
    for b=1:o.nBP
        UnbledBackground = zeros(o.nBP,o.nRounds);
        UnbledBackground(b,:) = 1;
        BackgroundEigenvectors(b,:,:) = o.z_scoreBleedMatrix(:,:,1) * UnbledBackground;
    end
    BackgroundEigenvectors = BackgroundEigenvectors./vecnorm(BackgroundEigenvectors(:,:),2,2);
    BackgroundEigenvalues = zeros(o.nBP,1);
else    
    if nargin<3 || isempty(DotProductThresh)
        DotProductThresh = 0.35;
    end
    if nargin<4 || isempty(nTilesToUse)
        nTilesToUse = 5;
    end
    NonemptyTiles = find(~o.EmptyTiles)';
    if size(NonemptyTiles,2)==1
        NonemptyTiles = NonemptyTiles';
    end
    fprintf('Finding background codes\n');
    %For speed, only use 5 tiles
    if isempty(o.BackgroundEigenvectorTiles)
        EigenvectorTiles = NonemptyTiles(randperm(numel(NonemptyTiles),min(nTilesToUse,length(NonemptyTiles))));
    else
        EigenvectorTiles = o.BackgroundEigenvectorTiles;
    end
    BackgroundSpotColors = cell(max(NonemptyTiles),1);
    
    for t=EigenvectorTiles
        [LocalYX,SpotColors] = o.get_spot_colors_all_pixels(t);
        SpotColors = (double(SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
        NormSpotColors = SpotColors./vecnorm(SpotColors(:,:),2,2);  %So norm=1        
        MaxGeneDotProduct = get_bled_code_max_dot_product(NormSpotColors,NormBledCodes,nCodes);
        KeepBackgroundSpots = MaxGeneDotProduct<DotProductThresh;
        BackgroundSpotColors{t} = SpotColors(KeepBackgroundSpots,:,:);
    end
    BackgroundSpotsFullSet = cell2mat(BackgroundSpotColors);
    CovMatrix = BackgroundSpotsFullSet(:,:)'*BackgroundSpotsFullSet(:,:);
    [Eig,BackgroundEigenvalues] = get_eig(CovMatrix);
    BackgroundEigenvectors = reshape(Eig',[o.nRounds*o.nBP,o.nBP,o.nRounds]);
    BackgroundEigenvalues = BackgroundEigenvalues/vecnorm(BackgroundEigenvalues);
end
%Find out how much background eigenvectors resemble genes
[BackgroundMaxGeneDotProduct,BackgroundMaxGeneDotProductGene] = ...
    get_bled_code_max_dot_product(BackgroundEigenvectors,NormBledCodes,nCodes);

if o.Graphics
    figure(180432); clf
    for i=1:min(9,length(BackgroundEigenvectors))
        subplot(3,3,i);
        imagesc(squeeze(BackgroundEigenvectors(i,:,:)));
        title(sprintf('Eigenvalue = %.3f\n%s DotProduct = %.3f',...
            BackgroundEigenvalues(i),o.GeneNames{BackgroundMaxGeneDotProductGene(i)},...
            BackgroundMaxGeneDotProduct(i)));
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'yTickLabel', o.bpLabels(1:o.nBP));
        set(gca, 'xtick', 1:o.nRounds);
        set(gca, 'XTickLabel', 1:o.nRounds);
        caxis([-1,1]);
        colormap(gca,bluewhitered);
        if i==4
            xlabel('Round')
            ylabel('Channel');
        end
    end
    drawnow;
    sgtitle('Background Eigenvectors');
end
     
end

function [MaxDotProduct,MaxGeneNumber] = get_bled_code_max_dot_product(SpotColors,GeneBledCodes,nCodes)
%For each spot, s, finds dot product of SpotColors(s,:) with each gene.
%Then returns max absolute dot product and corresponding gene. 
%Both SpotColors and GeneBledCodes should have been normalised. 
GeneDotProduct = zeros(size(SpotColors,1),nCodes);
for GeneNo=1:nCodes
    GeneDotProduct(:,GeneNo) = SpotColors(:,:)*GeneBledCodes(GeneNo,:)';
end
[MaxDotProduct,MaxGeneNumber] = max(abs(GeneDotProduct),[],2);
end   
