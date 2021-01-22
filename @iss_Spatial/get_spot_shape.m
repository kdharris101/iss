function o = get_spot_shape(o)
%% Find shape of spots to fit to data
% First, collect all isolated spots for all genes.
o.spShapeDiam = o.ExtractR2*3-1;
o.spShapeRad = (o.spShapeDiam-1)/2;
% have to detect peaks in image - all above thresh.
Small = 1e-6;
SE_find_peaks = strel('disk', o.PixelDetectRadius);

nCodes = length(o.CharCodes);
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end
AllGeneIsolatedIm = cell(nCodes,max(NonemptyTiles));     %Collect all genes together then get eigenvalues
fprintf('Finding spot shape\n');
for t=o.BackgroundEigenvectorTiles
    [LocalYX,SpotColors] = o.get_spot_colors_all_pixels(t);
    LocalYX = double(LocalYX);
    SpotColors = (double(SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    for GeneNo=1:nCodes
        %Turn coefficients into an image
        GeneIm = zeros(max(LocalYX));
        Ind = sub2ind(size(GeneIm),LocalYX(:,1),LocalYX(:,2));
        GeneIm(Ind) = SpotColors(:,:)*o.z_scoreBledCodes(GeneNo,:)';
        %Find peaks
        Dilate = imdilate(GeneIm, SE_find_peaks);
        MaxProjectPeakInd = find(GeneIm + Small >= Dilate & GeneIm > o.spShapeThresh);
        LocalInd = find(ismember(Ind,MaxProjectPeakInd));
        
        %find those that are isolated - far away from neighbour
        k0 = KDTreeSearcher(LocalYX(LocalInd,:));
        [~, d2] = k0.knnsearch(LocalYX(LocalInd,:), 'k', 2);    %find distance to neighbour
        SpotIsolated = d2(:,2)>o.spShapeIsolationDist;
        nIsolated = sum(SpotIsolated);
        IsolatedYX = LocalYX(LocalInd(SpotIsolated),:);
        AllIsolatedIm = zeros(nIsolated,o.spShapeDiam*o.spShapeDiam);
        for s=1:nIsolated
            YX = IsolatedYX(s,:);
            try
                SmallIm = GeneIm(YX(1)-o.spShapeRad:YX(1)+o.spShapeRad,YX(2)-o.spShapeRad:YX(2)+o.spShapeRad);
            catch
                %set to nan if overlaps with edge
                SmallIm = nan(o.spShapeDiam,o.spShapeDiam);
            end
            AllIsolatedIm(s,:) = SmallIm(:);
        end
        %get rid of nan
        AllIsolatedIm = AllIsolatedIm(~isnan(AllIsolatedIm(:,1)),:);
        AllGeneIsolatedIm{GeneNo,t} = AllIsolatedIm;
    end
end

%Form covariance matrix from spot images and get eigenvectors
%MAYBE MAKE DIFFERENT EIGENVECTORS FOR GENES AND BACKGROUNDS.
%PROBABLY UNNECCESSARY, IN WHICH CASE, DO THIS PART WHEN MAKE BACKGROUND
%EIGENVECTORS SO DONT HAVE TO LOAD PIXELS TWICE. EIGENVECTORS LOOK THE
%SAME. ALSO, MAYBE ONLY FIND THESE FROM SUBSET OF TILES.
AllGeneIsolatedIm = cell2mat(reshape(AllGeneIsolatedIm,1,[])');
%MAYBE NORM ALLGENEISOLATEDIM BEFORE TAKING COVMATRIX
CovMatrix = AllGeneIsolatedIm'*AllGeneIsolatedIm;
[o.spShapeEigenvectors,spShapeEigenvalues] = get_eig(CovMatrix);
o.spShapeEigenvalues = spShapeEigenvalues/vecnorm(spShapeEigenvalues);

if o.Graphics
    figure(180432); clf
    caxis_lim = [min(o.spShapeEigenvectors(:)),max(o.spShapeEigenvectors(:))];
    for i=1:min(9,size(o.spShapeEigenvectors,2))
        subplot(3,3,i);
        imagesc(reshape(o.spShapeEigenvectors(:,i),[o.spShapeDiam,o.spShapeDiam]));
        title(sprintf('Eigenvalue = %.3f',...
            o.spShapeEigenvalues(i)));
        set(gca, 'ytick', []);
        set(gca, 'yTickLabel', []);
        set(gca, 'xtick', []);
        set(gca, 'XTickLabel', []);
        caxis(caxis_lim);
        colormap(gca,bluewhitered);
    end
    sgtitle('Spot Shape Eigenvectors');
    drawnow;
end

end

