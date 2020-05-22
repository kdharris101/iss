function iss_view_spot(o, FigNo, ImSz, SpotLocation,ScoreMethod, SpotNum)
    %Check PCR by plotting location of spot in each round and color channel
    %FigNo is Figure number of plot3D if open otherwise set to any number
    %ImSz is radius of image that is plotted for each round and channel.
    %Default value is 7 pixels.
    %If SpotLocation is true, will use location of spot closest to
    %crosshair, otherwise will use actual position of crosshair. Default is
    %false.
    %ScoreMethod='DotProduct means use o.SpotCodeNo and ScoreMethod='Prob' means use
    %o.pSpotCodeNo to highlight Gene squares. Set to 'Prob' by default
    %SpotNum is index of spot that you want to look at.    
    
    %If Dist to nearest spot is less than MaxDist, then the crosshair will
    %be red on the squares corresponding to the gene the nearest spot
    %matches to.
    MaxDist = 10;
    
    if nargin<3 || isempty(ImSz)
        ImSz = 7;
    end
    if ImSz>100
        warning('ImSz too large, setting to 7');
        ImSz = 7;
    end
    
    if nargin<4 || isempty(SpotLocation)
        SpotLocation = false;
    end
    
    if nargin>=6
        SpotLocation = true;
        SpotNo = SpotNum;
        Dist = 0;
        if strcmpi('Pixel',ScoreMethod)
            xy = o.pxSpotGlobalYX(SpotNo,[2,1]);
        else
            xy = o.SpotGlobalYX(SpotNo,[2,1]);
        end
    else
        if nargin>=2
            figure(FigNo);
        end
        CrossHairColor = [1,1,1];   %Make white as black background
        xy = ginput_modified(1,CrossHairColor);
        S = evalin('base', 'issPlot2DObject');
        if nargin<5 || isempty(ScoreMethod)
            ScoreMethod = S.CallMethod;
        elseif ~strcmpi(S.CallMethod,ScoreMethod)
            if strcmpi('Prob',ScoreMethod) || strcmpi('DotProduct',ScoreMethod)
                S.SpotYX = o.SpotGlobalYX;
            elseif strcmpi('Pixel',ScoreMethod)
                S.SpotYX = o.pxSpotGlobalYX;
            end
            S.QualOK = 1;
        end
        InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
        PlotSpots = find(InRoi & S.QualOK);         %Only consider spots that can be seen in current plot
        [Dist,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
        SpotNo = PlotSpots(SpotIdx);
        if SpotLocation || round(Dist)==0
            SpotLocation = true;
            Dist = 0;
            xy = S.SpotYX(SpotNo,[2,1]);
        end

    end

    
    fprintf('loading channel/round images...');
    %Find tile that the point is on and local centered coordinates in reference round
    Dist2Tiles = [xy(2),xy(1)]-o.TileOrigin(:,:,o.ReferenceRound);
    Dist2Tile = min(sum(Dist2Tiles(Dist2Tiles(:,1)>=0 & Dist2Tiles(:,2)>=0,:),2));
    t = find(sum(Dist2Tiles,2)==Dist2Tile);
    LocalYX = [xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound)-o.TileCentre;
        
    if strcmpi(ScoreMethod,'DotProduct')
        numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.SpotCodeNo(SpotNo))),'\d','match'))+1;
    elseif strcmpi(ScoreMethod,'Prob')
        numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.pSpotCodeNo(SpotNo))),'\d','match'))+1;
    elseif strcmpi(ScoreMethod,'Pixel')
        numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.pxSpotCodeNo(SpotNo))),'\d','match'))+1;
    end

    
    try
        clf(27642)
        figure(27642)
    catch
        figure(27642)
    end
    set(gcf,'Position',[164,108,1621,805])
    for r=1:o.nRounds
        
        Ylegends = {o.bpLabels{:}};
        Xlegends = string(1:o.nRounds);
        for b=1:o.nBP
            
            rbYX = round(o.A(b)*([LocalYX,1]*o.D(:,:,t,r))+o.TileCentre);
            y0 = rbYX(1);
            x0 = rbYX(2);
            if y0>o.TileSz || y0<1 || x0>o.TileSz || x0<1
                continue;
            end
            y1 = max(1,y0 - ImSz);
            y2 = min(o.TileSz,y0 + ImSz);
            x1 = max(1,x0 - ImSz);
            x2 = min(o.TileSz,x0 + ImSz);
            BaseIm = int32(imread(o.TileFiles{r,t}, b, 'PixelRegion', {[y1 y2], [x1 x2]}))-o.TilePixelValueShift;
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
            axis([x0-ImSz, x0+ImSz, y0-ImSz, y0+ImSz]);
            if Dist<MaxDist
                caxis([0,max(150,max(BaseImSm(:)))]);
            else
                caxis([min(-50,min(BaseImSm(:))),max(150,max(BaseImSm(:)))]);
            end
            colorbar;
            if numCharCode(r)==b && Dist<MaxDist
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
    if SpotLocation
        if strcmpi(ScoreMethod,'DotProduct')
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
        elseif strcmpi(ScoreMethod,'Prob')
            %Color different parameters depending if over threshold
            if o.pSpotScore(SpotNo)>o.pScoreThresh
                c1 = [0,0.7,0]; else; c1 = [0,0,0];end
            if o.pLogProbOverBackground(SpotNo)<o.pLogProbThresh
                c2 = [1,0,0]; else; c2 = [0,0,0];end
            if o.pSpotScore(SpotNo)+o.pSpotScoreDev(SpotNo)<o.pDevThresh
                c3 = [1,0,0]; else; c3 = [0,0,0];end
            if o.pSpotIntensity(SpotNo)<o.pIntensityThresh
                c4 = [1,0,0]; else; c4 = [0,0,0];end
            figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
            figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProb = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
                SpotNo,o.GeneNames{o.pSpotCodeNo(SpotNo)},'\color[rgb]',c1,o.pSpotScore(SpotNo),'\color[rgb]',c2, o.pLogProbOverBackground(SpotNo),...
                '\color[rgb]',c3,o.pSpotScoreDev(SpotNo),'\color[rgb]',c4,o.pSpotIntensity(SpotNo));
        elseif strcmpi(ScoreMethod,'Pixel')
            %Color different parameters depending if over threshold
            if o.pxSpotScore(SpotNo)>o.pScoreThresh
                c1 = [0,0.7,0]; else; c1 = [0,0,0];end
            if o.pxLogProbOverBackground(SpotNo)<o.pLogProbThresh
                c2 = [1,0,0]; else; c2 = [0,0,0];end
            if o.pxSpotScore(SpotNo)+o.pxSpotScoreDev(SpotNo)<o.pDevThresh
                c3 = [1,0,0]; else; c3 = [0,0,0];end
            if o.pxSpotIntensity(SpotNo)<o.pIntensityThresh
                c4 = [1,0,0]; else; c4 = [0,0,0];end
            figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
            figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProb = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
                SpotNo,o.GeneNames{o.pxSpotCodeNo(SpotNo)},'\color[rgb]',c1,o.pxSpotScore(SpotNo),'\color[rgb]',c2, o.pxLogProbOverBackground(SpotNo),...
                '\color[rgb]',c3,o.pxSpotScoreDev(SpotNo),'\color[rgb]',c4,o.pxSpotIntensity(SpotNo));
        end
    end
end
    


