function iss_view_prob_loc(o,LookupTable,xy,ImSz)
%ISS_VIEW_PROB_LOC Summary of this function goes here
%   Detailed explanation goes here
xy = round(xy);
y_range = xy(2)-ImSz:xy(2)+ImSz;
x_range = xy(1)-ImSz:xy(1)+ImSz;
[A,B] = meshgrid(y_range,x_range);
c=cat(2,A',B');
GlobalYX=reshape(c,[],2);

t = o.get_local_tile(GlobalYX);
[RoundTile,~] = get_SpotTileEachRound(o,GlobalYX,t);
LocalYX = GlobalYX-o.TileOrigin(t,:,o.ReferenceRound);
load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
SpotColors = get_spot_colors(o,LocalYX,t,...
    RoundTile,AllBaseLocalYX);
AllLogProbOverBackground = o.get_LogProbOverBackground(SpotColors,LookupTable);
SortedLogProb = sort(AllLogProbOverBackground,2,'descend');
SecondBestLogProb = SortedLogProb(:,2);

g = 50;
ImYX = GlobalYX-min(GlobalYX)+1;
LogProbIm = zeros(max(ImYX));
ImInd = sub2ind(size(LogProbIm),ImYX(:,1),ImYX(:,2));
LogProbIm(ImInd) = AllLogProbOverBackground(:,g);
figure(57392);
imagesc(x_range,y_range,LogProbIm);
set(gca, 'YDir', 'normal');
colormap(gca,bluewhitered);


Small = 1e-6;
se1 = strel('disk', o.PixelDetectRadius);     %Needs to be bigger than in detect_spots
Dilate = imdilate(LogProbIm, se1);
MaxPixels = find(LogProbIm + Small >= Dilate);
PeakInd = find(ismember(ImInd,MaxPixels));
hold on
scatter(GlobalYX(PeakInd,2),GlobalYX(PeakInd,1),30,'gx','LineWidth',2);

InRoi = all(int64(round(o.pxSpotGlobalYX))>=min(GlobalYX) & round(o.pxSpotGlobalYX)<=max(GlobalYX),2);
CorrectGene = o.pxSpotCodeNo==g;
QualOK = quality_threshold(o,'Pixel');
scatter(o.pxSpotGlobalYX(InRoi&CorrectGene&QualOK,2),o.pxSpotGlobalYX(InRoi&CorrectGene&QualOK,1),...
    30,'yo','LineWidth',2);
scatter(o.pxSpotGlobalYX(InRoi&CorrectGene&~QualOK,2),o.pxSpotGlobalYX(InRoi&CorrectGene&~QualOK,1),...
    30,'yx','LineWidth',2);
hold off
end

