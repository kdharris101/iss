function [k, v, s2] = ScaledKMeans(x, v0)
% [k, v] = ScaledKMeans(x, v0);
%
% does a clustering that minimizes the norm of x_i - g_i * v_{k_i}, where:
% x_i is the i'th row of input matrix x (nPoints by nDims)
% g_i is a gain (not explicitly computed or returned)
% v_k is the k'th row of output v containing cluster means (nClusters by nDims)
% k_i is the i'th entry of output k giving the cluster of each point (nPoints by 1)
% s2_k is the first eigenval of the outer product matrix for cluster k
%
% input v0 (nClusters by nDims) is the starting point. (Required.)

ScoreThresh = 0; % only keep good matches 
MinClusterSize = 10; % delete clusters with too few points
ConvergenceCriterion = 0; % if this many or less changed, terminate

% normalize starting points and original data
vNorm = bsxfun(@rdivide, v0, sqrt(sum(v0.^2,2)));
v = vNorm;
xNorm = bsxfun(@rdivide, x, sqrt(sum(x.^2,2)));

[nClusters, nDims] = size(v);
s2 = zeros(nClusters,1);


MaxIter = 100;
k = nan; % to make sure it doesn't finish on first iteration

for i=1:MaxIter;
    kOld = k;
    
    score = xNorm*v'; % project each point onto each cluster. Use normalized so we can interpret score
    [TopScore, k] = max(score,[],2); % find best cluster for each point
    k(TopScore<ScoreThresh)=0; % unclusterable points
    
    if sum(k~=kOld)<=ConvergenceCriterion % need a better criterion!
        break;
    end

    % find top svd component for points assigned to each cluster
    for c=1:nClusters
        MyPoints = x(k==c,:); % don't use normalized, to avoid overweighting weak points
        nMyPoints = length(MyPoints);
        if nMyPoints<MinClusterSize
            v(c,:) = 0;
            continue;
        end
        [TopEvec, s2(c)] = eigs(double(MyPoints'*MyPoints)/nMyPoints, 1);
        v(c,:) = TopEvec*sign(mean(TopEvec)); % make them positive
    end
    

%     figure(9046);
%     hist(k,0:max(k));
%     drawnow

    
end

