function SpotNo = iss_view_codes(o, FigNo, Norm, SpotNum)
    %Now have option for different normalisations
    %Norm = 1: Raw Colors
    %Norm = 2: Normalised by percentile for each round and color channel
    %and then so each round has unit norm
    %Norm = 3: Normalised by percentile for each color channel across all
    %rounds

    
    if nargin>=4
        SpotNo = SpotNum;
    else
        if nargin>=2
            figure(FigNo);
        end
        CrossHairColor = [1,1,1];   %Make white as black background
        xy = ginput_modified(1,CrossHairColor);
        S = evalin('base', 'issPlot2DObject');
        if strcmpi('Pixel',S.CallMethod)
            S.SpotYX = o.SpotGlobalYX;
            S.QualOK = 1;
        end
        InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
        PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
        [~,SpotIdx] = min(sum(abs(o.SpotGlobalYX(PlotSpots,:)-[xy(2),xy(1)]),2));
        SpotNo = PlotSpots(SpotIdx);
    end
    CodeNo = o.SpotCodeNo(SpotNo);
    
    if nargin<3 || isempty(Norm)
        Norm = 1;
    end
    
    %Different Normalisations
    if isempty(Norm) || Norm == 1
        cSpotColors = o.cSpotColors;
        cBledCodes = o.pBledCodes;
        if isempty(o.pBledCodes)
            cBledCodes = o.BledCodes;
        end
    elseif Norm == 2
        cSpotColors = o.cNormSpotColors;
        cBledCodes = o.NormBledCodes;
    elseif Norm == 3
        cSpotColors = o.cSpotColors;
        NewBleedMatrix = o.pBleedMatrix;
        for b = 1:o.nBP
            bSpotColors = o.cSpotColors(:,b,:);
            p = prctile(bSpotColors(:), o.SpotNormPrctile);
            cSpotColors(:,b,:) = cSpotColors(:,b,:)/p;
            NewBleedMatrix(b,:,:) = o.pBleedMatrix(b,:,:)/p;                        
        end
        cBledCodes = change_bled_codes(o,NewBleedMatrix);
    end
    
    
    MeasuredCode = squeeze(cSpotColors(SpotNo,:,:));
    CodeShape = size(MeasuredCode);
    
    %Get square outlining unbled code
    gUnbled = reshape(o.UnbledCodes(CodeNo(1),:,:),CodeShape);
    gSquares = zeros(o.nRounds,4);
    for r=1:o.nRounds
        try
            gSquares(r,:) = [r-0.5,find(gUnbled(:,r,:)==1)-0.5,1,1];
        end
    end
    
    figure(930476530)
    subplot(2,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Measured code: match %.3f to %s', o.SpotScore(SpotNo), o.GeneNames{CodeNo}));
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    hold on
    for r=1:o.nRounds
        rectangle('Position',gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
    end
    hold off
    
    subplot(2,1,2)
    BledCode = cBledCodes(CodeNo,:);
    imagesc(reshape(BledCode, CodeShape)); colorbar
    caxis([0 max(BledCode(:))]);
    title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');
    hold on
    for r=1:o.nRounds
        rectangle('Position',gSquares(r,:),'EdgeColor','r','LineWidth',1,'LineStyle',':')
    end
    hold off

    fprintf('Spot %d at yx=(%d,%d): code %d, %s\n', ...
        SpotNo, o.SpotGlobalYX(SpotNo,1),o.SpotGlobalYX(SpotNo,2), CodeNo, o.GeneNames{CodeNo});

    
end
    