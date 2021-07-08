function o = PointCloudRegister6_GT(o, y0, x0, nTiles)     %MADE A THE SAME FOR ALL TILES
% o = o.PointCloudRegister(y, x, A0, Options)
% 
% Perform point cloud registration to map points x onto points y by
% iterative closest point: repeatedly finding the best y for each x, 
% and doing linear regression to find the M that maps best maps x to y
%
% inputs:
% y0 is a cell containig the YX location of all spots in all rounds 
% and colour channels for all tiles
%
% x0{t,b} is a cell containing the YX location of spots in the 
% reference round for tile t, channel b
%
% ToPlot: array of form [t,b,r] of specific example case to show plot of
% for debugging purposes
%
%%
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end
if isempty(o.TileCentre)
    o.TileCentre = 0.5*[o.TileSz+1,o.TileSz+1];
end
UseChannelsOrig = o.UseChannels;

%Centre SpotYX
x0(NonemptyTiles,o.ReferenceSpotChannels) = cellfun(@(x0) x0(:,1:2)-o.TileCentre,...
    x0(NonemptyTiles,o.ReferenceSpotChannels),'UniformOutput',false);
x = cell(nTiles,o.nBP);
for t=NonemptyTiles
    for b = o.ReferenceSpotChannels        
        %Append array of ones for translation
        x(t,b) = {[x0{t,b},ones(size(x0{t,b},1),1)]};
    end
end


if isempty(o.PcDist)
    o.PcDist = inf;
end

%Initialize variables
D = o.D;
for t=NonemptyTiles
    for r = o.UseRounds
        for b=o.UseChannels
            if r==o.AnchorRound&&(b==o.DapiChannel||b==o.AnchorChannel)
                D(:,:,t,r,b) = eye(3,2);
            else
                D(1:2,:,t,r,b) = eye(2);
                D(1,1,t,r,b) = o.A(b,1);
                D(2,2,t,r,b) = o.A(b,2);
                D(3,:,t,r,b) = o.D0(t,:,r);
            end
        end
    end
end

fprintf('\nPCR - Finding well isolated points');
% find well isolated points as those whose second neighbor is far
y = y0;
for t=NonemptyTiles
    for r=o.UseRounds
        for b=o.UseChannels
            
            % make kd tree - default options!
            k0 = KDTreeSearcher(y0{t,b,r});
            [~, d2] = k0.knnsearch(y0{t,b,r}, 'k', 2);
            if isfinite(o.PcDist) && size(y0{t,b,r},1) > 1 
                y(t,b,r) = {y0{t,b,r}(d2(:,2)>o.PcDist*2,:)};
            end
            
        end
    end
end

fprintf('\nPCR - Making kd trees');
%Make kd trees out of these well isolated points
k = cell(nTiles,o.nBP,o.nRounds);
for t=NonemptyTiles
    for r=o.UseRounds
        for b=o.UseChannels
            k(t,b,r) = {KDTreeSearcher(y{t,b,r})};
        end
    end
end


%%
UseMe = cell(nTiles,o.nBP,o.nRounds+o.nExtraRounds);           %nP DIFFERENT FOR DIFFERENT TILES!!!
Neighbor = cell(nTiles,o.nBP,o.nRounds+o.nExtraRounds);
MyNeighb = cell(nTiles,o.nBP,o.nRounds+o.nExtraRounds);
xM = cell(nTiles,o.nBP,o.nRounds+o.nExtraRounds);
nMatches = zeros(nTiles,o.nBP,o.nRounds+o.nExtraRounds);
Error = zeros(nTiles,o.nBP,o.nRounds+o.nExtraRounds);
IsConverged = false(nTiles,o.nBP,o.nRounds+o.nExtraRounds);
TotalNeighbMatches = length(NonemptyTiles)*length(o.UseChannels)*...
    length(o.UseRounds);
FinishedGoodImages = false;
i_finished_good = inf;

for i=1:o.PcIter+100
    
    LastNeighbor = Neighbor;
    
    
    vertcat(o.RawLocalYX{:,b});
    
    for t=NonemptyTiles
        x_t = vertcat(x{t,:});
        for r=o.UseRounds                
            for b=o.UseChannels
                xM(t,b,r) = {x_t*D(:,:,t,r,b)+o.TileCentre};
            end
        end
    end
           
    for t=NonemptyTiles
        x_t = vertcat(x{t,:});
        for r=o.UseRounds     
            if r==o.AnchorRound
                %If anchor forget dapi and anchor channels.
                o.UseChannels = setdiff(UseChannelsOrig,[o.AnchorChannel,o.DapiChannel]);
            else
                o.UseChannels = UseChannelsOrig;
            end
            for b=o.UseChannels
                if IsConverged(t,b,r); continue; end
                Neighbor(t,b,r) = {k{t,b,r}.knnsearch(xM{t,b,r})};
                [~,Dist] = k{t,b,r}.knnsearch(xM{t,b,r});
                UseMe(t,b,r) = {Dist<o.PcDist};
                MyNeighb(t,b,r) = {Neighbor{t,b,r}(UseMe{t,b,r}>0)};
                nMatches(t,b,r) = sum(UseMe{t,b,r});
                Error(t,b,r) = sqrt(mean(Dist(UseMe{t,b,r}>0).^2));
                if ~FinishedGoodImages
                    D(:,:,t,r,b) = x_t(UseMe{t,b,r}>0,:)\(y{t,b,r}(MyNeighb{t,b,r},:)-o.TileCentre);
                elseif FinishedGoodImages
                    %For transforms found that are outliers, find D again
                    %but now regularize against the median transform.
                    Scale = [o.PcRegLambda2;o.PcRegLambda2;o.PcRegLambda1]; %As shift differences are bigger than scaling differences.
                    %Need to equalise their influence.
                    x_regularised = [x_t(UseMe{t,b,r}>0,:);eye(3).*Scale];
                    y_regularised = [(y{t,b,r}(MyNeighb{t,b,r},:)-o.TileCentre); D_average(:,:,t,r,b).*Scale];
                    D(:,:,t,r,b) = x_regularised\y_regularised;
                end
            end
        end
    end
    o.UseChannels = UseChannelsOrig;
        
    if isempty(o.ToPlot) == 0
        t = o.ToPlot(1);
        b = o.ToPlot(2);
        r = o.ToPlot(3);
        if i == 1
            fprintf('\nPlotting tile %d, color channel %d, round %d', t, b,r);
        end
        figure(29387648);
        fprintf('\nIteration %d: %d matches, mean error %f', i, nMatches(t,b,r), Error(t,b,r));
        clf; hold on
        plot(y{t,b,r}(:,2), y{t,b,r}(:,1), 'g+');
        plot(xM{t,b,r}(:,2), xM{t,b,r}(:,1), 'r+');
        plot([xM{t,b,r}(UseMe{t,b,r}>0,2) y{t,b,r}(MyNeighb{t,b,r},2)]',...
            [xM{t,b,r}(UseMe{t,b,r}>0,1) y{t,b,r}(MyNeighb{t,b,r},1)]', 'k-', 'linewidth', 1);
        
        drawnow;
    end
    
    IsConverged = cellfun(@isequal, Neighbor, LastNeighbor);
    nNeighbMatches = sum(sum(sum(IsConverged(NonemptyTiles,o.UseChannels,o.UseRounds))));
    
    if (nNeighbMatches == TotalNeighbMatches && ~FinishedGoodImages) ||...
            (i==o.PcIter && ~FinishedGoodImages)
        [D, D_average, PcFailed] = o.get_bad_transforms_GT(NonemptyTiles,nMatches,D);
        IsConverged(PcFailed) = false;
        nNeighbMatches = sum(sum(sum(IsConverged(NonemptyTiles,o.UseChannels,o.UseRounds))));
        FinishedGoodImages = true;
        i_finished_good = i;        
    elseif i == i_finished_good+1
        IsConverged(PcFailed) = false;
        nNeighbMatches = sum(sum(sum(IsConverged(NonemptyTiles,o.UseChannels,o.UseRounds))));
    end
        
    fprintf('\nPCR - Iteration %d: Converged images = %d/%d',i,nNeighbMatches,TotalNeighbMatches);
    %if min(min(min(cellfun(@isequal, Neighbor, LastNeighbor)))) == 1; break; end
    if nNeighbMatches == TotalNeighbMatches; break; end
end

fprintf('\n');
if nNeighbMatches<o.PcCovergedImgFrac*TotalNeighbMatches
    warning('\nPCR - Less than %d%% of images have converged',o.PcCovergedImgFrac*100);
end

%%
o.D = D;
o.gtRegInfo.nMatches = nMatches;
o.gtRegInfo.Error = Error;
o.gtRegInfo.nPcCovergedImg = nNeighbMatches/TotalNeighbMatches;
o.gtRegInfo.PcFailed = PcFailed;

