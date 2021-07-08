function [MaxLogProb,MaxScore] = iss_view_prob_loc(o,LookupTable,xy,GeneNumbers,ImSz,AllBaseLocalYX)
%ISS_VIEW_PROB_LOC Summary of this function goes here
%   Detailed explanation goes here
y_range = xy(2)-ImSz:xy(2)+ImSz;
x_range = xy(1)-ImSz:xy(1)+ImSz;
[A,B] = meshgrid(y_range,x_range);
c=cat(2,A',B');
GlobalYX=reshape(c,[],2);

t = o.get_local_tile(GlobalYX);
[RoundTile,~] = get_SpotTileEachRound(o,GlobalYX,t);
LocalYX = GlobalYX-o.TileOrigin(t,:,o.ReferenceRound);
if nargin<6 || isempty(AllBaseLocalYX)
    load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
end
SpotColors = get_spot_colors(o,LocalYX,t,...
    RoundTile,AllBaseLocalYX);
AllLogProbOverBackground = o.get_LogProbOverBackground(SpotColors,LookupTable);
[MaxLogProb,MaxIdx] = max(AllLogProbOverBackground(:,GeneNumbers));
SortedLogProb = sort(AllLogProbOverBackground,2,'descend');
SecondBestLogProb = SortedLogProb(:,2);
MaxScore = MaxLogProb-SecondBestLogProb(MaxIdx);

if nargout==0
    %Add button to flip to spot score.
    
    GlobalYX = round(GlobalYX);
    nGenes = length(GeneNumbers);
    CaxisLims = [min(min(AllLogProbOverBackground(:,GeneNumbers))),...
        max(MaxLogProb)];
    Fig = figure(57393);
    set(Fig,'Position',[242,257,350+(nGenes-1)*350,408]);
    i=1;
    for g=GeneNumbers
        subplot(1,nGenes,i);
        ImYX = GlobalYX-min(GlobalYX)+1;
        LogProbIm = zeros(max(ImYX));
        ImInd = sub2ind(size(LogProbIm),ImYX(:,1),ImYX(:,2));
        LogProbIm(ImInd) = AllLogProbOverBackground(:,g);
        imagesc(x_range,y_range,LogProbIm);
        caxis(CaxisLims);
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
        title(sprintf('%s',o.GeneNames{g}));
        if i==nGenes
            colorbar;
        end
        i=i+1;
    end
end
end

