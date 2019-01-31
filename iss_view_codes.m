function debug_code(o, FigNo)


    if nargin>=2
        figure(FigNo);
    end
    
    set(gca, 'Color', [1 1 1]*.2);
     xy = ginput(1);
     set(gca, 'color', 'k');
%    x = 5832; y = 7936;
    [~,SpotNo] = min(sum(abs(o.SpotGlobalYX-[xy(2),xy(1)]),2));
    CodeNo = o.SpotCodeNo(SpotNo);
    
    MeasuredCode = squeeze(o.cSpotColors(SpotNo,:,:));
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
    