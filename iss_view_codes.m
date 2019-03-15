function debug_code(o, FigNo, ~)


    if nargin>=2
        figure(FigNo);
    end
    
    %If 3 or more arguments, use normed SpotColors that are actually used
    %to determine spot scores
    if nargin>=3
        SpotColors = bsxfun(@rdivide, o.cSpotColors, prctile(o.cSpotColors, o.SpotNormPrctile));
        FlatSpotColors = SpotColors(:,:);
        o.SpotIntensity = sqrt(sum(FlatSpotColors.^2,2));
        NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, o.SpotIntensity);
        cSpotColors = reshape(NormFlatSpotColors,size(o.cSpotColors));
    else
        cSpotColors = o.cSpotColors;
    end
    
    
    set(gca, 'Color', [1 1 1]*.2);
     xy = ginput(1);
     set(gca, 'color', 'k');
%    x = 5832; y = 7936;
    [~,SpotNo] = min(sum(abs(o.SpotGlobalYX-[xy(2),xy(1)]),2));
    CodeNo = o.SpotCodeNo(SpotNo);
    
    MeasuredCode = squeeze(cSpotColors(SpotNo,:,:));
    CodeShape = size(MeasuredCode);
    
    figure(930476530)
    subplot(2,1,1);
    imagesc(MeasuredCode); colorbar
    caxis([0 max(MeasuredCode(:))]);
    title(sprintf('Measured code: match %.3f to %s', o.SpotScore(SpotNo), o.GeneNames{CodeNo}));
    
    set(gca, 'ytick', 1:4);
    set(gca, 'YTickLabel', o.bpLabels);
    
    subplot(2,1,2)
    cBledCode = o.BledCodes(CodeNo,:);
    imagesc(reshape(cBledCode, CodeShape)); colorbar
    caxis([0 max(cBledCode(:))]);

    title(sprintf('Predicted Code for %s, code #%d', o.GeneNames{CodeNo}, CodeNo));
    
    
    set(gca, 'ytick', 1:4);
    set(gca, 'YTickLabel', o.bpLabels);

    fprintf('Spot %d at yx=(%d,%d): code %d, %s\n', ...
        SpotNo, o.SpotGlobalYX(SpotNo,1),o.SpotGlobalYX(SpotNo,2), CodeNo, o.GeneNames{CodeNo});

    
end
    