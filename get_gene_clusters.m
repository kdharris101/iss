function ClusterUseSpots = get_gene_clusters(o,r,k,UseSpots)
%This finds the index of spots that are clustered, such that they are with a
%distance r pixels of at least k other spots. 
%Inputs
%r - max distance in pixels for point to be within cluster
%k - number of points within distance r that is required to constitute a cluster
%UseSpots - Logical array of length, nSpots = length(o.pSpotScore),
%only spots s such UseSpots(s)==1 will be considered. 
%Defaults to the usual probability thresholding but say if you
%wanted to consider all spots then UseSpots = ones(nSpots,1);
%Output
%SpotNumbers - spot indexes of spots that satisfy SpotsToUse and that are
%in clusters

if nargin<2 || isempty(r)
    r = 18;
end
if nargin<3 || isempty(k)
    k = 20;
end
if nargin<4 || isempty(UseSpots)
    UseSpots = o.quality_threshold_prob;
end

Numbers = find(UseSpots);
GoodSpots = o.SpotGlobalYX(Numbers,:);

ClusterIdx = rangesearch(GoodSpots,GoodSpots,r);
ClusterSizes = cell2mat(cellfun(@(x) length(x),ClusterIdx,'UniformOutput',false));
GoodIdx = ClusterIdx(ClusterSizes>=k);
GoodIdx = cell2mat({cat(2, GoodIdx{:})});
SpotNumbers = Numbers(unique(GoodIdx));
ClusterUseSpots = UseSpots & ismember(1:length(UseSpots),SpotNumbers)';
end
