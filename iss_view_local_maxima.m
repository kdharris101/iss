function Dist2Maxima = iss_view_local_maxima(o, FigNo, Norm, Method, SpotNum)
%% SpotNo = iss_view_prob(o, FigNo, Norm, Method, SpotNum)
%
% This function lets you view the spot code, gene code
% and ln(Prob) for the best match for a chosen spot.
%
% o: iss object
% FigNo: figure number (default, current figure)
% Norm: normalization option (described below)
% Method: 'Prob' or 'Pixel' to consider gene assignments given
% by o.pSpotCodeNo or o.pxSpotCodeNo respectively.
% SpotNum: number of spot to analyze (default is to click)
% SpotNo: returns the number of the spot analyzed.
%
% Norm = 1: Raw Colors
% Norm = 2: Normalised by o.SpotNormPrctile in each colour channel and round,
% then if o.CallSpotsCodeNorm='WholeCode', normalise so whole code has L2 norm
% of 1 but if o.CallSpotsCodeNorm='Round', normalise so each round has L2 norm of 1.


%%
if nargin>=5
    SpotNo = SpotNum;
else
    if nargin>=2
        figure(FigNo);
    end
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    S = evalin('base', 'issPlot2DObject');
    if nargin<4 || isempty(Method)
        if strcmpi('DotProduct',S.CallMethod)
            Method = 'Prob';
        else
            Method = S.CallMethod;
        end
    elseif ~strcmpi(S.CallMethod,Method)
        if strcmpi('Prob',Method)
            S.SpotYX = o.SpotGlobalYX;
        elseif strcmpi('Pixel',Method)
            S.SpotYX = o.pxSpotGlobalYX;
        end
        S.QualOK = 1;
    end
    InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
    PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
    [~,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);
    
    
end

%Different parameters for different methods
if strcmpi('Prob',Method)
    CodeNo = o.pSpotCodeNo(SpotNo);
    SpotColor = o.cSpotColors(SpotNo,:,:);
    SpotScore = o.pSpotScore(SpotNo);
    LogProbOverBackground = o.pLogProbOverBackground(SpotNo);
    SpotScoreDev = o.pSpotScoreDev(SpotNo);
    SpotIntensity = o.pSpotIntensity(SpotNo);
    SpotGlobalYX = o.SpotGlobalYX(SpotNo,:);
elseif strcmpi('Pixel',Method)
    CodeNo = o.pxSpotCodeNo(SpotNo);
    SpotColor = o.pxSpotColors(SpotNo,:,:);
    SpotScore = o.pxSpotScore(SpotNo);
    LogProbOverBackground = o.pxLogProbOverBackground(SpotNo);
    SpotScoreDev = o.pxSpotScoreDev(SpotNo);
    SpotIntensity = o.pxSpotIntensity(SpotNo);
    SpotGlobalYX = o.pxSpotGlobalYX(SpotNo,:);
else
    error('Spot calling method not valid, should be Prob or Pixel');
end

if nargin<3 || isempty(Norm)
    Norm = 1;
end


Dist2Tiles = SpotGlobalYX-o.TileOrigin(:,:,o.ReferenceRound);
Dist2Tile = min(sum(Dist2Tiles(Dist2Tiles(:,1)>=0 & Dist2Tiles(:,2)>=0,:),2));
t = find(sum(Dist2Tiles,2)==Dist2Tile);
LocalYX = SpotGlobalYX-o.TileOrigin(t,:,o.ReferenceRound)-o.TileCentre;

%% Find nearest local maxima on each round/channel
% Record Spot intensity - Local maxima intensity and distance to local
% maxima.
Diff2Maxima = zeros(1,o.nBP,o.nRounds);
Dist2Maxima = zeros(o.nBP,o.nRounds);
ImSz = 8;
SpotBaseImPos = ImSz+1;
se1 = strel('disk', o.DetectionRadius);
Small = 1e-6; % just a small number, for computing local maxima: shouldn't matter what it is
for r=1:o.nRounds
    for b=1:o.nBP
        %Read in small image about correct coordinates
        rbYX = round([LocalYX,1]*o.D(:,:,t,r,b)+o.TileCentre);
        y0 = rbYX(1);
        x0 = rbYX(2);
        y1 = max(1,y0 - ImSz);
        y2 = min(o.TileSz,y0 + ImSz);
        x1 = max(1,x0 - ImSz);
        x2 = min(o.TileSz,x0 + ImSz);
        try
            BaseIm = int32(imread(o.TileFiles{r,t}, b, 'PixelRegion', {[y1 y2], [x1 x2]}))-o.TilePixelValueShift; 
        catch
            BaseIm = zeros(17,17);
        end
        Dilate = imdilate(BaseIm, se1);
        MaxPixels = find(BaseIm + Small >= Dilate & BaseIm>0);
        [yPeak, xPeak] = ind2sub(size(BaseIm), MaxPixels);
        dist=bsxfun(@hypot,yPeak-SpotBaseImPos,xPeak-SpotBaseImPos);
        ClosestIdx = find(dist==min(dist));
        if length(ClosestIdx)>1
            %Take local maxima with greatest intensity
            [~,Best] = max(BaseIm(MaxPixels(find(dist==min(dist)))));
            ClosestIdx = ClosestIdx(Best);
        end
        
        if isempty(dist) || y0<1 || y0>o.TileSz || x0<1 || x0>o.TileSz
            Diff2Maxima(1,b,r) = 0;
            Dist2Maxima(b,r) = 100;
        else
            Diff2Maxima(1,b,r) = double(BaseIm(SpotBaseImPos,SpotBaseImPos)) -...
                double(BaseIm(MaxPixels(ClosestIdx)));
            Dist2Maxima(b,r) = -dist(ClosestIdx);
        end
    end
end


% %% Different Normalisations
% if isempty(Norm) || Norm == 1
%     cSpotColor = Diff2Maxima;
%     cBledCodes = o.pBledCodes;
% elseif Norm == 2
%     if strcmpi(o.BleedMatrixType,'Separate')
%         p = prctile(o.cSpotColors, o.SpotNormPrctile);
%     elseif strcmpi(o.BleedMatrixType,'Single')
%         p = zeros(1,o.nBP,o.nRounds);
%         for b = 1:o.nBP
%             bSpotColors = o.cSpotColors(:,b,:);
%             p(:,b,:) = prctile(bSpotColors(:), o.SpotNormPrctile);
%         end
%     end
%     cSpotColor = Diff2Maxima./p;
%     %cSpotColor = cSpotColor/sqrt(sum(cSpotColor(:).^2));
%     %cSpotColor = o.cNormSpotColors(SpotNo,:,:);
%     cBledCodes = bsxfun(@rdivide, o.BledCodes, sqrt(sum(o.BledCodes.^2,2)));
%     %cBledCodes = o.NormBledCodes;
% end
% 
% %% Plot
% MeasuredCode = squeeze(cSpotColor);
% CodeShape = size(MeasuredCode);
% BledCode = cBledCodes(CodeNo,:);
% 
% %Get square outlining unbled code
% gUnbled = reshape(o.UnbledCodes(CodeNo(1),:,:),CodeShape);
% gSquares = zeros(o.nRounds,4);
% for r=1:o.nRounds
%     try
%         gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
%     end
% end
% 
% try
%     clf(430476533)
%     figure(430476533)
% catch
%     figure(430476533)
% end
% subplot(3,1,1);
% imagesc(MeasuredCode); colorbar
% %caxis([0 max(MeasuredCode(:))]);
% title(sprintf('Spot Code'));
% set(gca, 'ytick', 1:o.nBP);
% set(gca, 'YTickLabel', o.bpLabels);
% ylabel('Color Channel');
% hold on
% for r=1:o.nRounds
%     rectangle('Position',gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
% end
% hold off
% 
% subplot(3,1,2)
% imagesc(Dist2Maxima); colorbar
% %caxis([0 max(cBledCode(:))]);
% title('Distance to Local Maxima');
% set(gca, 'ytick', 1:o.nBP);
% set(gca, 'YTickLabel', o.bpLabels);
% ylabel('Color Channel');
% hold on
% for r=1:o.nRounds
%     rectangle('Position',gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
% end
% 
% subplot(3,1,3)
% imagesc(reshape(BledCode, CodeShape)); colorbar
% %caxis([0 max(cBledCode(:))]);
% title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
% set(gca, 'ytick', 1:o.nBP);
% set(gca, 'YTickLabel', o.bpLabels);
% ylabel('Color Channel');
% hold on
% for r=1:o.nRounds
%     rectangle('Position',gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
% end
% 
% %Color different parameters depending if over threshold
% if SpotScore>o.pScoreThresh
%     c1 = [0,0.7,0]; else; c1 = [0,0,0];end
% if LogProbOverBackground<o.pLogProbThresh
%     c2 = [1,0,0]; else; c2 = [0,0,0];end
% if SpotScore+SpotScoreDev<o.pDevThresh
%     c3 = [1,0,0]; else; c3 = [0,0,0];end
% if SpotIntensity<o.pIntensityThresh
%     c4 = [1,0,0]; else; c4 = [0,0,0];end
% 
% set(gcf,'Position',[350 100 1000 850])
% figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
% figtitle.String = sprintf('%s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProbOverBackground = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
%     '\color[rgb]',c1,SpotScore,'\color[rgb]',c2, LogProbOverBackground,...
%     '\color[rgb]',c3,SpotScoreDev,'\color[rgb]',c4,SpotIntensity);
% %figtitle.Color='red';
% drawnow
% 
% fprintf('Spot %d at yx=(%d,%d): code %d, %s\n', ...
%     SpotNo, SpotGlobalYX(1),SpotGlobalYX(2),...
%     CodeNo, o.GeneNames{CodeNo});
end
