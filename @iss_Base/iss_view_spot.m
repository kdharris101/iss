function iss_view_spot(o, FigNo, ImSz, SpotLocation, ScoreMethod, IncludeGT, Filter, SpotNum)
%% iss_view_spot(o, FigNo, ImSz, SpotLocation, ScoreMethod, IncludeGT, SpotNum)
%
% Check PCR by plotting location of spot in each round and color channel
%
% FigNo: figure number (default, current figure)
% ImSz: radius of image that is plotted for each round and channel.
% Default value is 7 pixels.
% SpotLocation: logical,  if true, will use location of spot closest to
% crosshair, otherwise will use actual position of crosshair. Default is false.
% ScoreMethod: 'DotProduct means use o.SpotCodeNo and ScoreMethod='Prob' means use
% IncludeGT: if true, will also plot the ground truth rounds.
% Filter: if true, will filter image (default), if false it will not and
% just return the focus stacked images.
% SpotCodeNo to highlight Gene squares. Set to 'Prob' by default
% SpotNum: index of spot that you want to look at.


%%
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

if nargin<6 || isempty(IncludeGT)
    IncludeGT = false;
end

if nargin<7 || isempty(Filter)
    Filter = true;
end

if nargin<4 || isempty(SpotLocation)
    SpotLocation = false;
end

if nargin>=8
    if length(SpotNum)==2
        SpotLocation = false;
        xy = [SpotNum(2),SpotNum(1)];
        S.SpotYX = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX']);
        [Dist,SpotNo] = min(sum(abs(S.SpotYX-[xy(2),xy(1)]),2));
    else
        SpotLocation = true;
        SpotNo = SpotNum;
        Dist = 0;
        xy = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX'])(SpotNo,[2,1]);
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
        S.SpotYX = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX']);
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

if ~ismember({ScoreMethod},o.CallMethods)
    error('Method invalid, must be member of o.CallMethods.');
end
pf = o.CallMethodPrefix(ScoreMethod);
numCharCode = str2double(regexp(cell2mat(o.CharCodes(o.([pf,'SpotCodeNo'])(SpotNo))),'\d','match'))+1;

fprintf('loading channel/round images...');
if SpotLocation == false
    %Find tile that the point is on and local centered coordinates in reference round
    t = o.get_local_tile([xy(2),xy(1)]);
else
    t = o.([pf,'LocalTile'])(SpotNo);
end
LocalYX = [xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound);

if IncludeGT
    o.UseRounds = [o.UseRounds,o.gtRounds];
    numCharCode = [numCharCode,990*o.gtRounds];
end
nRounds = max(o.UseRounds);
[RoundTile,~] = get_SpotTileEachRound(o,[xy(2),xy(1)],t);
load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
[SpotColor,PointCorrectedLocalYX] = get_spot_colors(o,LocalYX,t,...
    RoundTile,AllBaseLocalYX);
if SpotLocation==true
    if max(max(abs(double(o.([pf,'SpotColors'])(SpotNo,:,o.UseRounds(o.UseRounds<=o.nRounds)))...
            -SpotColor(:,:,o.UseRounds(o.UseRounds<=o.nRounds)))))>0
        warning('Spot Color found is different from than in o object');
    end
end

if ~Filter
    t_rawdata_round = zeros(max(o.UseRounds));
    for r=o.UseRounds
        t_rawdata_round(r) = str2double(o.TileFiles{1,RoundTile(r)}(end-4));
    end
end

try
    clf(27642)
    figure(27642)
catch
    figure(27642)
end
set(gcf,'Position',[164,108,1621,805]);
Ylegends = {o.bpLabels{:}};
Xlegends = string(1:o.nRounds);
if IncludeGT
    for i=1:length(o.gtRounds)
        Xlegends = [Xlegends,["gt"+string(i)]];
    end
end
for r=o.UseRounds  
    if ~Filter
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);  %raw data file name for round r
        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        bfreader.setId(imfile);
        % get some basic image metadata
        [nSeries, nSerieswPos, ~, nZstacks, ~, ~] = get_ome_tilepos(bfreader);
        scene = nSeries/nSerieswPos;
        bfreader.close();
        
        % a new reader per worker
        bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
        % use the memo file cached before
        bfreader.setId(imfile);
        bfreader.setSeries(scene*t_rawdata_round(r)-1);
    end
    for b=1:o.nBP
        
        rbYX = round(PointCorrectedLocalYX(1,:,r,b));
        y0 = rbYX(1);
        x0 = rbYX(2);
        if y0>o.TileSz || y0<1 || x0>o.TileSz || x0<1
            continue;
        end
        y1 = max(1,y0 - ImSz);
        y2 = min(o.TileSz,y0 + ImSz);
        x1 = max(1,x0 - ImSz);
        x2 = min(o.TileSz,x0 + ImSz);
        
        if Filter
            BaseIm = int32(imread(o.TileFiles{r,RoundTile(r)}, b, 'PixelRegion', {[y1 y2], [x1 x2]}))-o.TilePixelValueShift;
            if o.SmoothSize
                SE = fspecial3('ellipsoid',o.SmoothSize);
                BaseImSm = imfilter(BaseIm, SE);
            else
                BaseImSm = BaseIm;
            end
        else
            I = cell(nZstacks,1);
            for z = 1:nZstacks
                iPlane = bfreader.getIndex(z-1, b-1, 0)+1;
                I{z} = bfGetPlane(bfreader, iPlane, x1, y1, 2*ImSz+1, 2*ImSz+1);
                %I{z} = bfGetPlane(bfreader, iPlane);
            end
            % focus stacking
            BaseImSm = o.fstack_modified(I(o.FirstZPlane:end));
            %BaseImSm = BaseImSm(y1:y2,x1:x2);
        end
        
        h = subplot(o.nBP, nRounds, (b-1)*nRounds + r);
        if r == 1 && b == 1
            Pos1 = get(h,'position');
        elseif r == 1 && b == o.nBP
            Pos2 = get(h,'position');
        elseif r == nRounds && b == o.nBP
            Pos3 = get(h,'position');
        end
        imagesc([x1 x2], [y1 y2], BaseImSm); hold on
        %caxis([min(-150,min(BaseImSm(:))),max(150,max(BaseImSm(:)))]);
        if Filter
            colormap(gca,bluewhitered);
        else
            try
                caxis([min(BaseImSm(:)),5000]);
            catch
            end
        end
        axis([x0-ImSz, x0+ImSz, y0-ImSz, y0+ImSz]);
        if Dist<MaxDist
            %caxis([0,max(150,max(BaseImSm(:)))]);
        else
            %caxis([min(-50,min(BaseImSm(:))),max(150,max(BaseImSm(:)))]);
        end
        colorbar;
        if numCharCode(r)==b && Dist<MaxDist 
            ax = gca;
            ax.XColor = 'r';
            ax.YColor = 'r';
            plot(xlim, [y0 y0], 'g'); plot([x0 x0], ylim, 'g');
        elseif IncludeGT && o.gtGeneNo(r,b)>0
            ax = gca;
            if Dist<MaxDist && o.([pf,'SpotCodeNo'])(SpotNo)==o.gtGeneNo(r,b)
                ax.XColor = 'r';
                ax.YColor = 'r';
                plot(xlim, [y0 y0], 'g'); plot([x0 x0], ylim, 'g');
            else
                ax.XColor = 'y';
                ax.YColor = 'y';
                plot(xlim, [y0 y0], 'y'); plot([x0 x0], ylim, 'y');
            end
        else
            plot(xlim, [y0 y0], 'k'); plot([x0 x0], ylim, 'k');
        end
        
        if r==1; ylabel(Ylegends{b},'Color',[0.15 0.15 0.15]); end
        if b==o.nBP; xlabel(Xlegends(r),'Color',[0.15 0.15 0.15]); end
        %title(sprintf('Round %d, Base %d, Tile %d', r, b, t));
        %drawnow
        set(gca, 'YDir', 'normal');
        hold off
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
        if o.dpSpotScore(SpotNo)>o.CombiQualThresh
            c1 = [0,0.7,0]; else; c1 = [0,0,0];end
        if o.dpSpotScoreDev(SpotNo)<o.CombiDevThresh
            c2 = [1,0,0]; else; c2 = [0,0,0];end
        if o.dpSpotIntensity(SpotNo)<o.CombiIntensityThresh
            c3 = [1,0,0]; else; c3 = [0,0,0];end
        figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
        figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
            SpotNo,o.GeneNames{o.dpSpotCodeNo(SpotNo)},'\color[rgb]',c1,o.dpSpotScore(SpotNo),'\color[rgb]',c2, o.dpSpotScoreDev(SpotNo),...
            '\color[rgb]',c3,o.dpSpotIntensity(SpotNo));
    elseif strcmpi(ScoreMethod,'Prob') || strcmpi(ScoreMethod,'Pixel')
        %Color different parameters depending if over threshold
        if o.([pf,'SpotScore'])(SpotNo)>o.pScoreThresh
            c1 = [0,0.7,0]; else; c1 = [0,0,0];end
        if o.([pf,'LogProbOverBackground'])(SpotNo)<o.pLogProbThresh
            c2 = [1,0,0]; else; c2 = [0,0,0];end
        if o.([pf,'SpotScore'])(SpotNo)+o.([pf,'SpotScoreDev'])(SpotNo)<o.pDevThresh
            c3 = [1,0,0]; else; c3 = [0,0,0];end
        if o.([pf,'SpotIntensity'])(SpotNo)<o.pIntensityThresh
            c4 = [1,0,0]; else; c4 = [0,0,0];end
        figtitle = sgtitle('', 'interpreter', 'tex');   %'tex' required for colors
        figtitle.String = sprintf('Spot %.0f is %s: %s{%f %f %f}Score = %.1f, %s{%f %f %f}LogProb = %.0f, %s{%f %f %f}Score Deviation = %.1f, %s{%f %f %f}Intensity = %.0f',...
            SpotNo,o.GeneNames{o.([pf,'SpotCodeNo'])(SpotNo)},'\color[rgb]',c1,o.([pf,'SpotScore'])(SpotNo),'\color[rgb]',c2, o.([pf,'LogProbOverBackground'])(SpotNo),...
            '\color[rgb]',c3,o.([pf,'SpotScoreDev'])(SpotNo),'\color[rgb]',c4,o.([pf,'SpotIntensity'])(SpotNo));
    end
end
end



