function PCR_view(o,y0,x0,t_plot,b_plot,r_plot)
%% PCR_view(t,b,r)
%
% plots result of point cloud registration between (tile,channel,round) = 
% (t,o.ReferenceChannel,o.ReferenceRound) to (t,b,r)
%
% o: iss object
% y0 is a cell containig the YX location of all spots in all rounds 
% and colour channels for all tiles (AllBaseLocalYX in FindSpotsWorkspace.m)
% x0{t,b} is a cell containing the YX location of spots in the 
% reference round for tile t, channel b. (o.RawLocalYX)
% t_plot: tile to plot
% b_plot: channel
% r_plot: round
%% 
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';

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

% find well isolated points as those whose second neighbor is far
% only need to consider t,b,r of interest for plotting
y = y0;
k0 = KDTreeSearcher(y0{t_plot,b_plot,r_plot});
[~, d2] = k0.knnsearch(y0{t_plot,b_plot,r_plot}, 'k', 2);
if isfinite(o.PcDist) && size(y0{t_plot,b_plot,r_plot},1) > 1 
    y(t_plot,b_plot,r_plot) = {y0{t_plot,b_plot,r_plot}(d2(:,2)>o.PcDist*2,:)};
end

%Make kd trees out of these well isolated points
k = cell(nTiles,o.nBP,o.nRounds);
k(t_plot,b_plot,r_plot) = {KDTreeSearcher(y{t_plot,b_plot,r_plot})};

% Find neighbours and which ones are within the min required distance
Neighbor = cell(nTiles,o.nBP,o.nRounds);
UseMe = cell(nTiles,o.nBP,o.nRounds);      
MyNeighb = cell(nTiles,o.nBP,o.nRounds);
xM = cell(nTiles,o.nBP,o.nRounds);

x_t = vertcat(x{t_plot,:});
xM(t_plot,b_plot,r_plot) = {x_t*o.D(:,:,t_plot,r_plot,b_plot)+o.TileCentre}; 
Neighbor(t_plot,b_plot,r_plot) = {k{t_plot,b_plot,r_plot}.knnsearch(xM{t_plot,b_plot,r_plot})};
[~,Dist] = k{t_plot,b_plot,r_plot}.knnsearch(xM{t_plot,b_plot,r_plot});
UseMe(t_plot,b_plot,r_plot) = {Dist<o.PcDist};                
MyNeighb(t_plot,b_plot,r_plot) = {Neighbor{t_plot,b_plot,r_plot}(UseMe{t_plot,b_plot,r_plot}>0)};


figure(2987648);
clf;
%plot(y{t_plot,b_plot,r_plot}(:,2), y{t_plot,b_plot,r_plot}(:,1), 'g+');
%plot(xM{t_plot,b_plot,r_plot}(:,2), xM{t_plot,b_plot,r_plot}(:,1), 'r+');
anchor_NeighbYX = x_t(UseMe{t_plot,b_plot,r_plot}>0,1:2)+o.TileCentre;
anchor_transform_NeighbYX = xM{t_plot,b_plot,r_plot}(UseMe{t_plot,b_plot,r_plot}>0,:);     %Anchor neighbours
rb_NeighbYX = y{t_plot,b_plot,r_plot}(MyNeighb{t_plot,b_plot,r_plot},:);     %Imaging round neighbours
%plot([anchor_NeighbYX(:,2) rb_NeighbYX(:,2)]', [anchor_NeighbYX(:,1) rb_NeighbYX(:,1)]', 'k-', 'linewidth', 1);
quiver(rb_NeighbYX(:,2),rb_NeighbYX(:,1),anchor_transform_NeighbYX(:,2)-rb_NeighbYX(:,2),...
    anchor_transform_NeighbYX(:,1)-rb_NeighbYX(:,1),'color','k');
MaxShift = max(abs(squeeze(o.D0(t_plot,:,:))'));
xlim([1-MaxShift(2),o.TileSz+MaxShift(2)]);
ylim([1-MaxShift(1),o.TileSz+MaxShift(2)]);
title(sprintf('Tile %d, channel %d, round %d: %d matches, mean error %.3f',...
    t_plot,b_plot,r_plot, o.nMatches(t_plot,b_plot,r_plot), o.Error(t_plot,b_plot,r_plot)));
xlabel('X');
ylabel('Y');



x_shift = anchor_transform_NeighbYX(:,2)-rb_NeighbYX(:,2);
y_shift = anchor_transform_NeighbYX(:,1)-rb_NeighbYX(:,1);
edges = -o.PcDist:3:o.PcDist;

Sz = 100;
YMeanShift = zeros(o.TileSz-Sz,1);
XMeanShift = zeros(o.TileSz-Sz,1);
for i=1:o.TileSz-Sz
    ySpotsInRange = anchor_NeighbYX(:,1)>=i & anchor_NeighbYX(:,1)<i+Sz;
    xSpotsInRange = anchor_NeighbYX(:,2)>=i & anchor_NeighbYX(:,2)<i+Sz;
    YMeanShift(i) = mean(y_shift(ySpotsInRange));
    XMeanShift(i) = mean(x_shift(xSpotsInRange));
end

Coordinates = 1+Sz/2:o.TileSz-Sz/2;
figure(2741);
set(gcf,'Position',[74,108,1321,605]);
ax1 = subplot(3,3,[1,4,7]);
plot(Coordinates,YMeanShift);
hold on
plot(Coordinates,XMeanShift,'Color','#D95319');
hold off
ax1.XLabel.String = 'Relavent anchor round coordinate';
ax1.YLabel.String = 'Moving average error';
ax1.Title.String = sprintf('Y error as function of y coordiate\nX error as function of x coordinate');
legend('Y','X','Location','northwest');
xlim([0,o.TileSz]);

ax2 = subplot(3,3,2);
histogram(y_shift,edges,'Normalization','Probability');
ax2.Title.String = sprintf('Y error between neighbours with dist<%d\nMean = %.3f',o.PcDist,mean(y_shift));
ax3 = subplot(3,3,3);
histogram(x_shift,edges,'Normalization','Probability','FaceColor','#D95319');
ax3.Title.String = sprintf('X error between neighbours with dist<%d\nMean = %.3f',o.PcDist,mean(x_shift));
linkaxes([ax2, ax3], 'x');
linkaxes([ax2, ax3], 'y');


yshift_fit = polyfit(anchor_NeighbYX(:,1),y_shift,1);
xshift_fit = polyfit(anchor_NeighbYX(:,2),x_shift,1);
ax4 = subplot(3,3,[5,8]);
scatter(anchor_NeighbYX(:,1),y_shift,11);
hold on
plot(1:o.TileSz,yshift_fit(1)*(1:o.TileSz)+yshift_fit(2),'k');
hold off
xlim([0,o.TileSz]);
ax4.XLabel.String = 'Anchor Round Y';
ax4.YLabel.String = 'Y error';
ax5 = subplot(3,3,[6,9]);
scatter(anchor_NeighbYX(:,2),x_shift,11,'MarkerEdgeColor','#D95319');
hold on
plot(1:o.TileSz,xshift_fit(1)*(1:o.TileSz)+xshift_fit(2),'k');
hold off
xlim([0,o.TileSz]);
ax5.XLabel.String = 'Anchor Round X';
ax5.YLabel.String = 'X error';
linkaxes([ax1, ax4, ax5], 'x');


figtitle = sgtitle('');   %'tex' required for colors
figtitle.String = sprintf(sprintf('Tile %d, channel %d, round %d: %d matches, mean error %.3f',...
    t_plot,b_plot,r_plot, o.nMatches(t_plot,b_plot,r_plot), o.Error(t_plot,b_plot,r_plot)));


% % Considering variation in both dimensions
% Sz=350;
% YMeanShift = zeros(o.TileSz-Sz,o.TileSz-Sz);
% XMeanShift = zeros(o.TileSz-Sz,o.TileSz-Sz);
% for x=1:o.TileSz-Sz
%     for y=1:o.TileSz-Sz
%         SpotsInRange = rb_NeighbYX(:,2)>=x & rb_NeighbYX(:,2)<x+Sz &...
%         rb_NeighbYX(:,1)>=y & rb_NeighbYX(:,1)<y+Sz;
%         YMeanShift(y,x) = mean(y_shift(SpotsInRange));
%         XMeanShift(y,x) = mean(x_shift(SpotsInRange));
%     end
% end
% Coordinates = 1+Sz/2:o.TileSz-Sz/2;
% figure(2742);
% h_y = imagesc(YMeanShift);
% set(h_y, 'XData', Coordinates);
% set(h_y, 'YData', Coordinates);
% xlim([1+Sz/2,o.TileSz-Sz/2]);
% ylim([1+Sz/2,o.TileSz-Sz/2]);
% xlabel('rb X coordinate');
% ylabel('rb Y coordinate');
% title('Y mean error');
% colorbar;
% figure(2743);
% h_x = imagesc(XMeanShift);
% set(h_x, 'XData', Coordinates);
% set(h_x, 'YData', Coordinates);
% xlim([1+Sz/2,o.TileSz-Sz/2]);
% ylim([1+Sz/2,o.TileSz-Sz/2]);
% xlabel('rb X coordinate');lp
% ylabel('rb Y coordinate');
% title('X mean error');
% colorbar;

drawnow;
end