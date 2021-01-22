function iss_view_spatial(o, FigNo, ImSzView, SpotLocation, SpotNum)
%% iss_view_spatial(o, FigNo, ImSzView, SpotLocation, SpotNum)
%
% Check call_spots_spatial by plotting dot product images and convolution
% images for each gene centered on chosen location.
%
% FigNo: figure number (default, current figure)
% ImSz: radius of image that is plotted for each gene.
% Default value is 7 pixels.
% SpotLocation: logical,  if true, will use location of spot closest to
% crosshair, otherwise will use actual position of crosshair. Default is false.
% SpotNum: index of spot that you want to look at. Or yx global coordinate
% of interest.


%%
pf = o.CallMethodPrefix('Spatial'); 

if nargin<3 || isempty(ImSzView)
    ImSzView = 7;
end
if ImSzView>100
    warning('ImSz too large, setting to 7');
    ImSzView = 7;
end
ImSzFull = ImSzView+o.spShapeRad;     %So get valid convolution in ImSzView.

if nargin<4 || isempty(SpotLocation)
    SpotLocation = false;
end

if nargin>=5
    if max(size(SpotNum))==1
        SpotLocation = true;
        SpotNo = SpotNum;
        xy = o.pxSpotGlobalYX(SpotNo,[2,1]);
    else
        SpotLocation = false;
        xy = [SpotNum(2),SpotNum(1)];
        [Dist,SpotNo] = min(sum(abs(o.([pf,'SpotGlobalYX'])-SpotNum),2));
        if round(Dist)==0
            SpotLocation = true;
        end
    end
else
    if nargin>=2
        figure(FigNo);
    end
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    S_plot = evalin('base', 'issPlot2DObject');
    InRoi = all(int64(round(S_plot.SpotYX))>=S_plot.Roi([3 1]) & round(S_plot.SpotYX)<=S_plot.Roi([4 2]),2);
    PlotSpots = find(InRoi & S_plot.QualOK);         %Only consider spots that can be seen in current plot
    [Dist,SpotIdx] = min(sum(abs(S_plot.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);
    if SpotLocation || round(Dist)==0
        SpotLocation = true;
        xy = S_plot.SpotYX(SpotNo,[2,1]);
    end    
end

%Find tile that the point is on
if ~SpotLocation
    t = o.get_local_tile([xy(2),xy(1)]);
else
    t = o.([pf,'LocalTile'])(SpotNo);
end

%Get YX positions of pixels needed for tile t.
S = struct;     %All info needed for plotting saved in struct.
S.LocalY0X0 = round([xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound));     
%View range is local tile coordinates that are viewed at the end.
%min/max stuff is for dealing situation on edge of tile.
S.ViewRange.X = max(S.LocalY0X0(2)-ImSzView,1):min(S.LocalY0X0(2)+ImSzView,o.TileSz);
S.ViewRange.Y = max(S.LocalY0X0(1)-ImSzView,1):min(S.LocalY0X0(1)+ImSzView,o.TileSz);
%Below coordinates is the range of local pixel location on tile that
%consider when doing convolutions. Needs to be larger that above as spots
%outside view range can influence pixels in view range.
ExpandRangeFactor = 5;      %ExpandRangeFactor = 3 gave different Coefs between S.Peak and S.Close
LocalY_range = S.LocalY0X0(1)-ImSzFull*ExpandRangeFactor:S.LocalY0X0(1)+ImSzFull*ExpandRangeFactor;
LocalX_range = S.LocalY0X0(2)-ImSzFull*ExpandRangeFactor:S.LocalY0X0(2)+ImSzFull*ExpandRangeFactor;
[A,B] = meshgrid(LocalX_range,LocalY_range);
c=cat(2,A',B');
AnchorLocalXY = reshape(c,[],2);
AnchorLocalYX = AnchorLocalXY(:,[2,1]);   %So ascending in Y to match call_spots_spatial.
%get all spot colors required. 
[LocalYX,SpotColors] = o.get_spot_colors_all_pixels(t,AnchorLocalYX);
SpotColors = (double(SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;

%Save required plotting info to struct
%S.SpotColors = cell(o.sp_nIter+1,1);     %Save spot colors for each iteration.
%S.SpotColors{1} = SpotColors;       
S.LocalYX = LocalYX;
S.nCodes = size(o.spBledCodes,1);
S.nGeneCodes = length(o.CharCodes);
S.iter = 1;         %Background images are for before this iteration.
                    %Scatter plot shows spots found in this iteration.
S.Convolve = 0;     %Initial image is dot product.
S.BledCodes = o.spBledCodes(:,:);
%S.DotProduct(:,g) will be the initial plot for gene g. 
S.DotProduct = SpotColors(:,:)*S.BledCodes';      
S.GeneNames = o.GeneNames;
if S.nCodes>length(o.CharCodes)
    for g=length(o.CharCodes)+1:S.nCodes
        S.GeneNames{g}='Bckgrnd';
    end
end 

%Index of spot of interest in gene subplots i.e. where crosshair is.
%Subtelty when near edge.
S.x0 = find(S.ViewRange.X==S.LocalY0X0(:,2));       
S.y0 = find(S.ViewRange.Y==S.LocalY0X0(:,1));

%Collect convolution images for each iteration as well as the spots found
%in each iteration - these are saved as S.PeakSpots.
[AllSpotCodeNo,AllSpotCoef,AllSpotLocalInd,AllSpotIter,AllResidualReduction,...
    S.SpotColors,S.ConvImages] = o.get_spatial_spots(SpotColors,LocalYX);
AllLocalYX = LocalYX(AllSpotLocalInd,:);

%Only keep those close to centre
MaxDist = ceil(sqrt(ImSzFull^2+ImSzFull^2));
[ViewSpots,~] = rangesearch(double(AllLocalYX),S.LocalY0X0,MaxDist);
ViewSpots = cell2mat(ViewSpots);
S.PeakSpots.LocalYX = AllLocalYX(ViewSpots,:);
S.PeakSpots.Iter = AllSpotIter(ViewSpots,:);
S.PeakSpots.CodeNo = AllSpotCodeNo(ViewSpots);
S.PeakSpots.Coef = AllSpotCoef(ViewSpots,:);
S.PeakSpots.Score = AllResidualReduction(ViewSpots,:);
[~,SortedIdx]=sortrows([S.PeakSpots.Iter,S.PeakSpots.CodeNo]);  %So ordering same as S.CloseSpots
S.PeakSpots.LocalYX = S.PeakSpots.LocalYX(SortedIdx,:);
S.PeakSpots.Iter = S.PeakSpots.Iter(SortedIdx);
S.PeakSpots.CodeNo = S.PeakSpots.CodeNo(SortedIdx);
S.PeakSpots.Coef = S.PeakSpots.Coef(SortedIdx,:);
S.PeakSpots.Score = S.PeakSpots.Score(SortedIdx,:);

%Get closest spots found in pipeline (i.e. saved in o object). Then see if
%they match S.PeakSpots. If don't remove duplicates, they should be exactly
%the same. S.PeakSpots should include duplicates whereas S.CloseSpots
%should not hence size(S.PeakSpots) should be >= size(S.CloseSpots) and all
%S.CloseSpots should appear in S.PeakSpots.
S.CloseSpots = get_spots_in_range(o,xy,MaxDist,t,pf);
if size(S.CloseSpots.Coef,1) ~= size(S.PeakSpots.Coef,1)
    nSpotsFoundNow = size(S.PeakSpots.Coef,1);
    nSpotsFoundPipeline = size(S.CloseSpots.Coef,1);
    warning(sprintf('Now found %d spots, whereas in pipeline, found %d spots',...
        nSpotsFoundNow,nSpotsFoundPipeline));
elseif sum(sum(abs(S.CloseSpots.Coef-S.PeakSpots.Coef)>0.1))>0
    nWrong = sum(sum(abs(S.CloseSpots.Coef-S.PeakSpots.Coef)>0.1,2)>0);
    nTotal = size(S.CloseSpots.CodeNo,1);
    warning(sprintf('Coeficients in convolution are different now as to the pipeline for %d/%d spots',...
        nWrong,nTotal));   
end
    
try
    clf(276465)
    figure(276465);
catch
    figure(276465);
end
set(gcf,'Position',[164,108,1621,805]);
%Keep colorbar limits constant for all iterations.
S.DotProductClim = [min(S.DotProduct(:)),max(S.DotProduct(:))];
%S.ConvolveClim = [-log10(S.SpotConvThresh), ...
%    max(log10(S.ConvImages{1}(:)))-log10(S.SpotConvThresh)];
S.ConvolveClim = [min(S.ConvImages{1}(:)),max(S.ConvImages{1}(:))];
plot_gene_data(S);


S.button = uicontrol('style','togglebutton',...
                 'unit','pix',...
                 'position',[20 40 60 20],...
                 'val',S.Convolve,...
                 'String','Convolve',...
                 'Tag','button1','callback',{@convolve_image});
             
S.slider = uicontrol('style','slide',...
 'unit','pix',...
 'position',[100 40 1293 18],...
 'min',1,'max',o.sp_nIter,'val',S.iter,...
 'sliderstep',[1/(o.sp_nIter-1) 1/(o.sp_nIter-1)],...
 'callback',{@change_iter});
%Need to assign S to workspace so can exchange data between uicontrol
%functions.
assignin('base','issViewSpatialObject',S);
end

function CloseSpots = get_spots_in_range(o,xy,Sz,t,pf)
    %Finds all spots within Sz of xy in tile t.
    %First consider genes.
    nShapeUse = min(o.nShapeEigUse,size(o.([pf,'Coefs']),2));
    [Idx,~] = rangesearch(o.([pf,'SpotGlobalYX']),xy([2,1]),Sz);
    Idx = cell2mat(Idx);
    CloseSpots.LocalYX = o.([pf,'SpotGlobalYX'])(Idx,:)-o.TileOrigin(t,:,o.ReferenceRound);
    CloseSpots.Iter = o.([pf,'SpotIter'])(Idx);
    CloseSpots.CodeNo = o.([pf,'SpotCodeNo'])(Idx);
    CloseSpots.Coef = o.([pf,'Coefs'])(Idx,1:nShapeUse);
    CloseSpots.Score = o.([pf,'SpotScore'])(Idx);
    %Sort by iteration found and then by code number.
    [~,SortedIdx]=sortrows([CloseSpots.Iter,CloseSpots.CodeNo]);
    CloseSpots.LocalYX = CloseSpots.LocalYX(SortedIdx,:);
    CloseSpots.Iter = CloseSpots.Iter(SortedIdx);
    CloseSpots.CodeNo = CloseSpots.CodeNo(SortedIdx);
    CloseSpots.Coef = CloseSpots.Coef(SortedIdx,:);
    CloseSpots.Score = CloseSpots.Score(SortedIdx);
end

function plot_gene_data(S)
%Function to plot either dot product or convolution data for each gene.
NextBestGene = S.PeakSpots.CodeNo(S.PeakSpots.Iter==S.iter);
InNextIter = S.PeakSpots.Iter==S.iter;
PastIter = S.PeakSpots.Iter<S.iter;
FutureIter = S.PeakSpots.Iter>S.iter;
%Get spot coordinates in small subplot reference frame.
ScatterPlotYX = double(S.PeakSpots.LocalYX)-double(S.LocalY0X0)+[S.y0,S.x0];
if S.Convolve
    S.climits = S.ConvolveClim;
else
    S.climits = S.DotProductClim;
end
gray = [0.6,0.6,0.6];   %Color for all genes not found   
dark_green = [0 0.5 0]; %Color for genes found on previous iteration
Small = 1e-20;            %So don't take log(0);
for g=1:S.nCodes   
    subplot(10, 8, g);
    CodeIm = zeros(max(S.LocalYX));
    Ind = sub2ind(size(CodeIm),S.LocalYX(:,1),S.LocalYX(:,2));
    CodeIm(Ind) = S.DotProduct(:,g);     
    if S.Convolve
        %S.ConvImages{S.iter}(S.ConvImages{S.iter}<0)=0;     %HACK FOR NEGATIVES BUT NEGATIVES ARE TRUE WITH ONE EIGENVECTOR
        %DataIm = log10(S.ConvImages{S.iter}(:,:,g)+Small)-log10(S.SpotConvThresh);
        if g>S.nGeneCodes
            %No Conv images for background as do by pixel instead.  
            DataIm = S.ConvImages{S.iter}(:,:,1)*0;
        else
            DataIm = S.ConvImages{S.iter}(:,:,g);
        end
    else
        DataIm = CodeIm;
    end
    DataIm = DataIm(S.ViewRange.Y,S.ViewRange.X);
    imagesc(DataIm); hold on;
    caxis(S.climits);
    colormap(gca,bluewhitered);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
    ax = gca;
    %Different colour schemes for when spot found or not.
    if ismember(g,NextBestGene)
        title([num2str(g), ': ',S.GeneNames{g}],'FontWeight','normal','Color',dark_green);
        ax.XColor = dark_green;
        ax.YColor = dark_green;
        plot(xlim, [S.y0 S.y0], 'Color', dark_green); plot([S.x0 S.x0], ylim, 'Color', dark_green);        
    elseif ismember(g,S.PeakSpots.CodeNo)
        title([num2str(g), ': ',S.GeneNames{g}],'FontWeight','normal','Color','k');
        ax.XColor = 'k';
        ax.YColor = 'k';
        plot(xlim, [S.y0 S.y0], 'k'); plot([S.x0 S.x0], ylim, 'k');
    else
        title([num2str(g), ': ',S.GeneNames{g}],'FontWeight','normal','Color',gray);
        ax.XColor = gray;
        ax.YColor = gray;
        plot(xlim, [S.y0 S.y0], 'Color',gray); plot([S.x0 S.x0], ylim,'Color',gray);
    end
    scatter(ScatterPlotYX(S.PeakSpots.CodeNo==g&InNextIter,2),...
        ScatterPlotYX(S.PeakSpots.CodeNo==g&InNextIter,1),60,'gx');
    scatter(ScatterPlotYX(S.PeakSpots.CodeNo==g&FutureIter,2),...
        ScatterPlotYX(S.PeakSpots.CodeNo==g&FutureIter,1),60,'kx');
    scatter(ScatterPlotYX(S.PeakSpots.CodeNo==g&PastIter,2),...
        ScatterPlotYX(S.PeakSpots.CodeNo==g&PastIter,1),60,dark_green,'x');
    set(gca, 'YDir', 'normal');
    hold off
end
%Single colorbar for all subplots.
cbh = colorbar;
cbh.Position(3) = cbh.Position(3)*2;
cbh.Position(4) = cbh.Position(4)*10;
cbh.Position(1) = .95-cbh.Position(3)/2;
cbh.Position(2) = 0.5-cbh.Position(4)/2;
sgtitle(sprintf('Iteration %d: %d spots (Not all visible)',...
    S.iter,sum(InNextIter)));
end

function S=convolve_image(h,~)
%uicontrol function to convolve.
S = evalin('base', 'issViewSpatialObject');
S.Convolve = round(get(h,'value'));
assignin('base','issViewSpatialObject',S);
plot_gene_data(S);
end

function S=change_iter(h,~)
%uicontrol function to change iteration.
S = evalin('base', 'issViewSpatialObject');
S.iter = round(get(h,'value'));
S.DotProduct = S.SpotColors{S.iter}(:,:)*S.BledCodes';
assignin('base','issViewSpatialObject',S);
plot_gene_data(S);
end
