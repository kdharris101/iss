function o = PCR_analysis(o,y0,x0)
%% o = PCR_analysis(o,y0,x0)
%
% Fits straight line fit to distribution of PCR shifts with position. I.e.
% if chromatic aberration is a problem, we would expect a strong gradient.
%
% o: iss object
% y0 is a cell containig the YX location of all spots in all rounds 
% and colour channels for all tiles (AllBaseLocalYX in FindSpotsWorkspace.m)
% x0{t,b} is a cell containing the YX location of spots in the 
% reference round for tile t, channel b. (o.RawLocalYX)
% saves:
% o.PcGrad
% o.PcMean
%% 
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';

if ~isfield(o.FindSpotsInfo,'nMatchesFailed')
    o.FindSpotsInfo.nMatchesFailed = o.nMatches;
    o.FindSpotsInfo.ErrorFailed = o.Error;
end

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

y = y0;
k = cell(nTiles,o.nBP,o.nRounds);
% Find neighbours and which ones are within the min required distance
Neighbor = cell(nTiles,o.nBP,o.nRounds);
UseMe = cell(nTiles,o.nBP,o.nRounds);      
MyNeighb = cell(nTiles,o.nBP,o.nRounds);
xM = cell(nTiles,o.nBP,o.nRounds);


fprintf('PCR - Doing analysis of tile   ');
for t=NonemptyTiles
    x_t = vertcat(x{t,:});
    for r=o.UseRounds
        for b=o.UseChannels   
            
            if t<10
                fprintf('\b%d',t);
            else
                fprintf('\b\b%d',t);
            end
            
            % find well isolated points as those whose second neighbor is far
            % make kd tree - default options!
            k0 = KDTreeSearcher(y0{t,b,r});
            [~, d2] = k0.knnsearch(y0{t,b,r}, 'k', 2);
            if isfinite(o.PcDist) && size(y0{t,b,r},1) > 1
                y(t,b,r) = {y0{t,b,r}(d2(:,2)>o.PcDist*2,:)};
            end
            %Make kd trees out of these well isolated points
            k(t,b,r) = {KDTreeSearcher(y{t,b,r})};
            
            xM(t,b,r) = {x_t*o.D(:,:,t,r,b)+o.TileCentre};
            Neighbor(t,b,r) = {k{t,b,r}.knnsearch(xM{t,b,r})};
            [~,Dist] = k{t,b,r}.knnsearch(xM{t,b,r});
            UseMe(t,b,r) = {Dist<o.PcDist};
            MyNeighb(t,b,r) = {Neighbor{t,b,r}(UseMe{t,b,r}>0)};            
            %o.nMatches(t,b,r) = sum(UseMe{t,b,r});
            
            anchor_NeighbYX = x_t(UseMe{t,b,r}>0,1:2);      %position in anchor round
            anchor_transform_NeighbYX = xM{t,b,r}(UseMe{t,b,r}>0,:);     %transformed anchor neighbours
            rb_NeighbYX = y{t,b,r}(MyNeighb{t,b,r},:);     %Imaging round neighbours  
            y_shift = anchor_transform_NeighbYX(:,1)-rb_NeighbYX(:,1);
            x_shift = anchor_transform_NeighbYX(:,2)-rb_NeighbYX(:,2);
            
            %Do straight line fit to both Y and Xdata and save gradient and intercept.
            yshift_fit = polyfit(anchor_NeighbYX(:,1)/o.TileCentre(1),y_shift,1);
            xshift_fit = polyfit(anchor_NeighbYX(:,2)/o.TileCentre(2),x_shift,1);
            o.PcGrad(t,b,r,1) = yshift_fit(1);
            o.PcMean(t,b,r,1) = mean(y_shift);
            o.PcGrad(t,b,r,2) = xshift_fit(1);
            o.PcMean(t,b,r,2) = mean(x_shift);
            
        end
    end
end
fprintf('\n');


end