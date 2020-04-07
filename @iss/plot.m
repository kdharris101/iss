function plot(o, BackgroundImageFile, Roi, Method)
% o.plot(BackgroundImage, Roi)
%
% plot the results of in situ sequencing spot detection. 
% SpotYX gives coordinates; Gene is a list of strings; 
% 
% if BackgroundImage is empty or missing, loaded from file o.BigDapiFile
% and subsetted to Roi (if present)
% If it is a numerical array, that is plotted (not subsetted for ROI)
% If zero, nothing is plotted
%
% Roi = [xmin xmax ymin ymax] shows only this part. Whole thing
% shown if empty or missing. Must be integers, xmin and ymin must be 1
%
% If Method == 'Prob', this changes the gene assignments to those given by
% the probability method rather than the dot product. In this case
% o.pScoreThresh is the threshold.
%
% sizes can be a vector or a scalar - only used for scatter, which isn't
% called anyway.
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

if nargin<4 || isempty(Method)
    Method = 'DotProduct';
end

if nargin<3 || isempty(Roi)
    Roi = round([1, max(o.SpotGlobalYX(:,2)), ...
    1, max(o.SpotGlobalYX(:,1))]);
end


if Roi(1) ~= 1 || Roi(3) ~= 1
    %Causes bugs if doesn't start at 1
    warning('Set min Roi to 1');
    Roi(1) = 1;
    Roi(3) = 1;
end

if (nargin<2 || isempty(BackgroundImageFile)) && ~isempty(o.BigDapiFile)
    if exist(o.BigDapiFile, 'file')
        fprintf('loading background image...');
        BackgroundImageFile = o.BigDapiFile;
        %Load in Dapi image
        Image = imread(BackgroundImageFile,'PixelRegion', {Roi(3:4), Roi(1:2)});
        fprintf('done\n');
    else
        warning('not sure what to do with BackgroundImage, setting to off');
        Image = zeros(Roi(4),Roi(2),'uint16');
    end
    
elseif ~isempty(BackgroundImageFile) && ~isnumeric(BackgroundImageFile)
    if exist(BackgroundImageFile, 'file')
        %Load in Dapi image
        fprintf('loading background image...');
        Image = imread(BackgroundImageFile,'PixelRegion', {Roi(3:4), Roi(1:2)});
        fprintf('done\n');
     else
        warning('not sure what to do with BackgroundImage, setting to off');
        Image = zeros(Roi(4),Roi(2),'uint16');
    end
    
elseif isnumeric(BackgroundImageFile)
    try
        Image = BackgroundImageFile(Roi(3):Roi(4),Roi(1):Roi(2));
    catch
        warning('Background Image wrong size - setting to off');
        Image = zeros(Roi(4),Roi(2),'uint16');
    end    
end


try
    S.FigNo = 234321;
    clf(S.FigNo)
catch
    S.FigNo = 234321;
end
S.fh = figure(S.FigNo);set(S.fh,'units','pixels','position',[500 200 800 600]);  %Left, Bottom, Width, Height
set(gcf, 'color', 'k');
set(gca, 'color', 'k');

S.Image = Image;
S.Background = imagesc(S.Image); hold on; colormap bone;
%set(S.Background, 'XData', [Roi(3), Roi(4)]);
%set(S.Background, 'YData', [Roi(1), Roi(2)]);
xlim([Roi(1) Roi(2)]);
ylim([Roi(3) Roi(4)]);

%title(['Z Plane ' num2str(S.MinZ)],'Color','w');


set(gca, 'YDir', 'normal');
%axis on

if strcmpi('Prob',Method)
    S.SpotGeneName = o.GeneNames(o.pSpotCodeNo);
    S.uGenes = unique(S.SpotGeneName);
    S.SpotYX = o.SpotGlobalYX;
    S.CallMethod = 'Prob';
elseif strcmpi('Pixel',Method)
    S.SpotGeneName = o.GeneNames(o.pxSpotCodeNo);
    S.uGenes = unique(S.SpotGeneName);
    S.SpotYX = o.pxSpotGlobalYX;
    S.CallMethod = 'Pixel';
else
    S.SpotGeneName = o.GeneNames(o.SpotCodeNo);
    S.uGenes = unique(S.SpotGeneName);
    S.SpotYX = o.SpotGlobalYX;
    S.CallMethod = 'DotProduct';
end
S.QualOK = o.quality_threshold(S.CallMethod);

S.Roi = Roi(1:4);
InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));

%hold on; GET RID OF HOLD AND JUST DELETE PLOTS USING DELETE MEANS THAT THE
%ZOOM IS KEPT BETWEEN Z PLANES
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYX(MySpots,2), S.SpotYX(MySpots,1), '.');
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

assignin('base','issPlot2DObject',S)

S.sl = uicontrol('Style','text','callback',{@sl_call,S},'BackgroundColor',[0,0,0]);  
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )





function [] = sl_call(varargin)
% Callback for the slider.
S = evalin('base', 'issPlot2DObject');  %Take S from workspace so can use iss_change_plot
h = findobj('type','line'); %KEY LINES: DELETE EXISTING SCATTER PLOTS SO CHANGE_SYMBOLS WORKS
delete(h)

set(gca, 'YDir', 'normal');
%axis on
%Stitle(['Z Plane ' num2str(ZPlane)],'Color','w');
InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
PlotSpots = find(InRoi & S.QualOK);
[~, S.GeneNo] = ismember(S.SpotGeneName(PlotSpots), S.uGenes);
S.h = zeros(size(S.uGenes));
%hold on;
for i=1:length(S.uGenes)
    MySpots = PlotSpots(S.GeneNo==i);
    if any(MySpots)
        S.h(i) = plot(S.SpotYX(MySpots,2), S.SpotYX(MySpots,1), '.');
    end
end 
%hold off

legend(S.h(S.h~=0), S.uGenes(S.h~=0));
legend off;

set(gca, 'Clipping', 'off');

if ~isempty(PlotSpots)
    change_gene_symbols(0);
else
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');    
end

%Update current Z position in woprkspace so can use for iss_view_codes
assignin('base','issPlot2DObject',S);








