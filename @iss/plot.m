function plot(o, BackgroundImage, Roi)
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
% shown if empty or missing. Must be integers
%
% sizes can be a vector or a scalar - only used for scatter, which isn't
% called anyway.
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html


BackgroundImageLabel = 'Background_DAPI';

%% load background image and subset it
% first see if it is already in figure to save time loading
if (nargin<2 || isempty(BackgroundImage)) && ~isempty(o.BigDapiFile)
    fprintf('loading background image...');
    if nargin<3 || isempty(Roi)
        Dapi = imread(o.BigDapiFile);
        dSize = size(Dapi);
        Roi = [1 dSize(2) 1 dSize(1)];
    else
        Dapi = imread(o.BigDapiFile, 'PixelRegion', {Roi(3:4), Roi(1:2)});
    end
        
    fprintf('done\n');
elseif nargin>= 2 && isnumeric(BackgroundImage) 
    %Dapi = BackgroundImage;
    Dapi = BackgroundImage;
else
    warning('not sure what to do with BackgroundImage, setting to off');
    Dapi = 0;
end

if nargin<3 || isempty(Roi)
    if numel(Dapi)>1
        Roi = [1, size(Dapi,2), 1, size(Dapi,1)];
    else
        Roi = [min(o.SpotGlobalYX(:,2)), max(o.SpotGlobalYX(:,2)), ...
            min(o.SpotGlobalYX(:,1)), max(o.SpotGlobalYX(:,1))];
    end
end
        



%% now plot it
clf; set(gcf, 'color', 'k');
set(gca, 'color', 'k');

if numel(Dapi)>1
    if ~isempty(Roi)
        hDapi = imagesc(Roi(1:2), Roi(3:4), Dapi);
    end
    colormap bone;
    set(hDapi, 'UserData', BackgroundImageLabel);
end

hold on;
set(gca, 'YDir', 'normal');
axis on

%% now find which spots to plot

SpotGeneName = o.GeneNames(o.SpotCodeNo);
uGenes = unique(SpotGeneName);

% which ones pass quality threshold (combi first)
QualOK = o.quality_threshold;

% now show only those in Roi
if ~isempty(Roi)
    InRoi = all(o.SpotGlobalYX>=Roi([3 1]) & o.SpotGlobalYX<=Roi([4 2]),2);
else
    InRoi = true;
end

PlotSpots = find(InRoi & QualOK);
[~, GeneNo] = ismember(SpotGeneName(PlotSpots), uGenes);
h = zeros(size(uGenes));
for i=1:length(uGenes)
    MySpots = PlotSpots(GeneNo==i);
    if any(MySpots)
        h(i) = plot(o.SpotGlobalYX(MySpots,2), o.SpotGlobalYX(MySpots,1), '.');
    end
end 

% add and remove legend to add DisplayNames
legend(h(h~=0), uGenes(h~=0));
legend off;

set(gca, 'Clipping', 'off');

change_gene_symbols(0);

end

