function iss_view_spot(o, FigNo, ScoreMethod, SpotNum)
    %Check PCR by plotting location of spot in each round and color channel
    %FigNo is Figure number of plot3D if open otherwise set to any number
    %ScoreMethod=1 means use o.SpotCodeNo and ScoreMethod=2 means use
    %o.pSpotCodeNo to highlight Gene squares. Set to 2 by default
    %SpotNum is index of spot that you want to look at.    

    if nargin>=4
        SpotNo = SpotNum;
    else
        if nargin>=2
            figure(FigNo);
        end
        CrossHairColor = [1,1,1];   %Make white as black background
        xy = ginput_modified(1,CrossHairColor);
        [~,SpotNo] = min(sum(abs(o.SpotGlobalYX(:,1:2)-[xy(2),xy(1)]),2));
    end
    if nargin < 3
        ScoreMethod = 2;
    end
    
    %ndRoundTile(s,r) tells us the tile spot s in the anchor round appears in
    %round r.
    %ndPointCorrectedLocalYXZ(s,:,r,b) tells us the YXZ position of spot s
    %in round r, channel b on tile ndRoundTile(s,r).
    %Good tells us which spots have non NaN cSpotColors in every
    %round/channel
    fprintf('loading channel/round images...');
    load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'),'ndRoundTile','ndPointCorrectedLocalYX','Good');
    GoodRoundTile = ndRoundTile(Good,:);
    GoodCorrectedYXZ = ndPointCorrectedLocalYX(Good,:,:,:);
    [nY, nX] = size(o.EmptyTiles);
    plsz = 7;
        
    if ScoreMethod == 1 
        numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.SpotCodeNo(SpotNo))),'\d','match'))+1;
    else
        numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.pSpotCodeNo(SpotNo))),'\d','match'))+1;
    end

    
    try
        clf(27642)
        figure(27642)
    catch
        figure(27642)
    end
    set(gcf,'Position',[164,108,1621,805])
    for r=1:o.nRounds
        t=GoodRoundTile(SpotNo,r);
        
        Ylegends = {o.bpLabels{:}};
        Xlegends = string(1:o.nRounds);
        for b=1:o.nBP
            
            y0 = GoodCorrectedYXZ(SpotNo,1,r,b);
            x0 = GoodCorrectedYXZ(SpotNo,2,r,b);
            if ~isfinite(x0) || ~isfinite(y0)
                continue;
            end
            y1 = max(1,y0 - plsz);
            y2 = min(o.TileSz,y0 + plsz);
            x1 = max(1,x0 - plsz);
            x2 = min(o.TileSz,x0 + plsz);
            BaseIm = int16(imread(o.TileFiles{r,t}, b, 'PixelRegion', {[y1 y2], [x1 x2]}))-o.TilePixelValueShift;
            if o.SmoothSize
                SE = fspecial3('ellipsoid',o.SmoothSize);
                BaseImSm = imfilter(BaseIm, SE);
            else
                BaseImSm = BaseIm;
            end
            
            h = subplot(o.nBP, o.nRounds, (b-1)*o.nRounds + r);
            if r == 1 && b == 1
                Pos1 = get(h,'position');
            elseif r == 1 && b == o.nBP
                Pos2 = get(h,'position');
            elseif r == o.nRounds && b == o.nBP
                Pos3 = get(h,'position');
            end
            imagesc([x1 x2], [y1 y2], BaseImSm); hold on
            axis([x0-plsz, x0+plsz, y0-plsz, y0+plsz]);
            caxis([0,max(150,o.cSpotColors(SpotNo,b,r))]);
            colorbar;
            if numCharCode(r)==b
                ax = gca;
                ax.XColor = 'r';
                ax.YColor = 'r';
                plot(xlim, [y0 y0], 'r'); plot([x0 x0], ylim, 'r');
            else
                plot(xlim, [y0 y0], 'w'); plot([x0 x0], ylim, 'w');
            end
            
            if r==1; ylabel(Ylegends{b},'Color',[0.15 0.15 0.15]); end
            if b==o.nBP; xlabel(Xlegends(r),'Color',[0.15 0.15 0.15]); end
            %title(sprintf('Round %d, Base %d, Tile %d', r, b, t));
            %drawnow
        end
    end
    PosDev = 0.02;
    SuperAxisPos = [Pos2(1:2)-PosDev,Pos3(1)+Pos3(2)-Pos2(1)+PosDev*2,Pos1(2)+Pos1(4)-Pos3(2)+PosDev*2];
    hSuper=axes('position',SuperAxisPos,'visible','off'); 
    hSuper.XLabel.Visible='on';
    hSuper.YLabel.Visible='on';
    axes(hSuper);
    ylabel('Channel');
    xlabel('Round');
    
    fprintf('done\n');
    
    if ScoreMethod == 1
        if o.SpotScore(SpotNo)>o.CombiQualThresh
            c1 = [0,0.7,0]; else; c1 = [0,0,0];end
        if o.SpotScoreDev(SpotNo)<o.CombiDevThresh
            c2 = [1,0,0]; else; c2 = [0,0,0];end
        if o.SpotIntensity(SpotNo)<o.CombiIntensityThresh
            c3 = [1,0,0]; else; c3 = [0,0,0];end
        figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
        figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
            SpotNo,o.GeneNames{o.SpotCodeNo(SpotNo)},'\color[rgb]',c1,o.SpotScore(SpotNo),'\color[rgb]',c2, o.SpotScoreDev(SpotNo),...
            '\color[rgb]',c3,o.SpotIntensity(SpotNo));        
    else
        %Color different parameters depending if over threshold
        if o.pSpotScore(SpotNo)>o.pScoreThresh
            c1 = [0,0.7,0]; else; c1 = [0,0,0];end
        if o.pLogProb(SpotNo)<o.pLogProbThresh
            c2 = [1,0,0]; else; c2 = [0,0,0];end
        if o.pSpotScore(SpotNo)+o.pSpotScoreDev(SpotNo)<o.pDevThresh
            c3 = [1,0,0]; else; c3 = [0,0,0];end
        if o.pSpotIntensity(SpotNo)<o.pIntensityThresh
            c4 = [1,0,0]; else; c4 = [0,0,0];end        
        figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
        figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProb = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
            SpotNo,o.GeneNames{o.pSpotCodeNo(SpotNo)},'\color[rgb]',c1,o.pSpotScore(SpotNo),'\color[rgb]',c2, o.pLogProb(SpotNo),...
            '\color[rgb]',c3,o.pSpotScoreDev(SpotNo),'\color[rgb]',c4,o.pSpotIntensity(SpotNo));
    end
    
end
    


