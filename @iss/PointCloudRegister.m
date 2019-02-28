function o = PointCloudRegister_NoAnchor(o, y, x0, A0, t, Options)
% [M, Error] = PointCloudRegister(y, x, M0, DistScale, Options)
% 
% Perform point cloud registration to map points x onto points y by
% iterative closest point: repeatedly finding the best y for each x, 
% and doing linear regression to find the M that maps best maps x to y
%
% inputs:
% y is a cell containig the centered YX location of all spots in all rounds 
% and colour channels
%
% x0 is the non centered YX location of spots in the anchor channel
%
% A0 are the initial scaling matrices for each colour channel 
% taking account of chromatic aberration. All default to identity if not
% specified
% 
% D0 are the initial shifts between rounds. All default to [0,0] if not
% specified
%
% DistScale: any x whose nearest neighbor is further than this won't count,
% default inf
%
% Options: what type of fit. For now ignored, the only option is a general linear
% model where x gets an extra column of ones and M is 2x3.
%%
MaxIter = 100;
Interactive = 0; % 1 to show, 2 to pause

[nP, nD] = size(x0);

%Centre x coordinates
x = minus(x0,[o.TileSz/2,o.TileSz/2]);

if nargin<4 || isempty(A0)
    A0 = zeros(nD,nD,o.nBP);
    for b=1:o.nBP
        A0(:,:,b) = eye(nD);
    end
elseif max(size(A0))==nD
    A = zeros(nD,nD,o.nBP);
    for b=1:o.nBP
        A(:,:,b) = A0;
    end
    A0 = A;    
end
A = A0;


if isempty(o.D0(:,:,t))
    o.D0(:,:,t) = zeros(o.nRounds,2);
end
D = o.D0(:,:,t);


if isempty(o.PcDist)
    o.PcDist = inf;
end


for r=1:o.nRounds 
    for b=1:o.nBP
        % make kd tree - default options!
        k0 = KDTreeSearcher(y{r,b});

        % find well isolated points as those whose second neighbor is far
        [~, d2] = k0.knnsearch(y{r,b}, 'k', 2);
        if isfinite(o.PcDist) && size(y{r,b},1) > 1 
            y(r,b) = {y{r,b}(d2(:,2)>o.PcDist*2,:)};
        end
    end
end

k = cell(o.nRounds,o.nBP);
for r=1:o.nRounds 
    for b=1:o.nBP
        k(r,b) = {KDTreeSearcher(y{r,b})};
    end
end

%%
Neighbor = zeros(nP,o.nRounds,o.nBP);
UseMe = zeros(nP,o.nRounds,o.nBP);
MyNeighb = cell(o.nRounds,o.nBP);
xM = zeros(nP,2,o.nRounds,o.nBP);
nMatches = zeros(o.nRounds,o.nBP);
Error = zeros(o.nRounds,o.nBP);

for i=1:MaxIter
    
    LastNeighbor = Neighbor;
    for r=1:o.nRounds 
        for b=1:o.nBP
            xM(:,:,r,b) = (A(:,:,b)*plus(x,D(r,:))')';   
        end
    end
    
    %This part finds new neighbours and new estimates for A
    for b=1:o.nBP
        xA = [];
        yA = [];
        for r=1:o.nRounds           
            [Neighbor(:,r,b), Dist] = k{r,b}.knnsearch(xM(:,:,r,b));
            UseMe(:,r,b) = Dist<o.PcDist;
            nMatches(r,b) = sum(UseMe(:,r,b));
            MyNeighb(r,b) = {Neighbor(UseMe(:,r,b)>0,r,b)};
            Error(r,b) = sqrt(mean(Dist(UseMe(:,r,b)>0).^2));
            
            xShift = plus(x(UseMe(:,r,b)>0,:),D(r,:));         %Add shift between rounds here
            xA = vertcat(xA, xShift);      
            yA = vertcat(yA, y{r,b}(MyNeighb{r,b},:));                     
            
        end
        A(:,:,b) = xA\yA;
    end
    
    %This part finds new estimates for D
    for r=1:o.nRounds
        xD = [];
        yD = [];
        for b=1:o.nBP
            xD = vertcat(xD,x(UseMe(:,r,b)>0,:));
            yScaled = (inv(A(:,:,b))*y{r,b}(MyNeighb{r,b},:)')';
            yD = vertcat(yD, yScaled);
        end
        D(r,:) = mean(minus(yD,xD));
    end
    
            
   
%    if Interactive
%        figure(29387648);
%        fprintf('Iteration %d: %d matches, mean error %f\n', i, sum(UseMe), Error);
%        clf; hold on
%         plot(y(:,2), y(:,1), 'g+');
%         plot(xM(:,2), xM(:,1), 'r+');
%        plot([xM(UseMe,2) y(MyNeighb,2)]', [xM(UseMe,1) y(MyNeighb,1)]', 'w-', 'linewidth', 1);
%
%        drawnow;
%        if Interactive>=2
%             pause
%        end
%    end
    
    if isequal(LastNeighbor, Neighbor); break; end
    
end

%%
o.A(:,:,:,t) = A;
o.D(:,:,t) = D;
o.nMatches(:,:,t) = nMatches;
o.Error(:,:,t) = Error;

