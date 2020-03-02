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
        CrossHairColor = [1,1,1];       %Make white as black background
        xy = ginput_modified(1,CrossHairColor);
        [~,SpotNo] = min(sum(abs(o.SpotGlobalYX-[xy(2),xy(1)]),2));
    end
    CodeNo = o.SpotCodeNo(SpotNo);
    
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
    
    figure(930476530)
    subplot(2,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Measured code: match %.3f to %s', o.SpotScore(SpotNo), o.GeneNames{CodeNo}));
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    
    subplot(2,1,2)
    BledCode = cBledCodes(CodeNo,:);
    imagesc(reshape(BledCode, CodeShape)); colorbar
    caxis([0 max(BledCode(:))]);

    title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
    
    
    set(gca, 'ytick', 1:o.nBP);
    set(gca, 'YTickLabel', o.bpLabels);
    ylabel('Color Channel');
    xlabel('Round');

    fprintf('Spot %d at yxz=(%d,%d): code %d, %s\n', ...
        SpotNo, o.SpotGlobalYX(SpotNo,1),o.SpotGlobalYX(SpotNo,2), CodeNo, o.GeneNames{CodeNo});

    
end
    