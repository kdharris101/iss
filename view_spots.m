%Check PCR by plotting spots found on anchor channel image and
%corresponding transformed points on specified base image.

function view_spots(o, tile, r, ColorChannel, ~)

    [nY, nX] = size(o.EmptyTiles);
    t = tile;
    b = ColorChannel;       %start at 1
    rr = o.ReferenceRound;
    
    %Get anchor round spots - only ones that can be read for all rounds
    AnchorYX = bsxfun(@minus, o.SpotGlobalYX, o.TileOrigin(t,:,rr));
    InTile = all(AnchorYX'>= 0 & AnchorYX' <=o.TileSz )';
    TileYX = AnchorYX(InTile,:);
    
    %Load anchor round image
    [y,x] = ind2sub([nY nX], t);
    AnchorIm = imread(o.TileFiles{rr,y,x}, o.AnchorChannel);
    if o.SmoothSize
        AnchorImSm = imfilter(AnchorIm, fspecial('disk', o.SmoothSize));
    else
        AnchorImSm = AnchorIm;
    end
    
    %Plot anchor round image and spots
    figure(72750);
    clf; set(gcf, 'color', 'k');
    imagesc(AnchorImSm);
    colormap bone;
    hold on
    set(gca, 'YDir', 'normal');
    plot(TileYX(:,2), TileYX(:,1), 'r.');
    hold off
    
    %Get position of spots on round r, colour channel b
    CenteredTileYX = TileYX - [o.TileSz/2,o.TileSz/2];
    CenteredTileCorrectedYX = (o.A(:,:,b)*(CenteredTileYX + o.D(t,:,r))')';
    AllMyPointCorrectedYX = round(CenteredTileCorrectedYX + [o.TileSz/2,o.TileSz/2]);
    InTile = all(AllMyPointCorrectedYX'>= 0 & AllMyPointCorrectedYX' <=o.TileSz )';
    MyPointCorrectedYX = AllMyPointCorrectedYX(InTile,:);
    
    %load image of round r, colour channel b
    FileName = o.TileFiles{r,t};
    TifObj = Tiff(FileName);
    TifObj.setDirectory(o.FirstBaseChannel + b - 1);
    BaseIm = TifObj.read();
    if o.SmoothSize
        BaseImSm = imfilter(double(BaseIm), fspecial('disk', o.SmoothSize));
    else
        BaseImSm = BaseIm;
    end
    
    % if give extra input, plots points without shift - can compare more
    % easily
    if nargin>=5
        figure(72751);
        clf; set(gcf, 'color', 'k');
        imagesc(BaseImSm);
        colormap bone;
        hold on
        set(gca, 'YDir', 'normal');
        plot(TileYX(:,2), TileYX(:,1), 'r.');
        hold off
    end
    
    %plot shifted points against BaseIm
    figure(72752);
    clf; set(gcf, 'color', 'k');
    imagesc(BaseImSm);
    colormap bone;
    hold on
    set(gca, 'YDir', 'normal');
    plot(MyPointCorrectedYX(:,2), MyPointCorrectedYX(:,1), 'r.');
    hold off
              
end
    


