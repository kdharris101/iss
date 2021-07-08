function [] = viewfilter_plot()
%% [] = viewfilter_plot()
%
% Once raw images for tile t and round r, loaded by view_filtering(o,r,t).
% This plots the results. Initial images are raw data, then you can press
% filter button to show filtered image. You can also change the radius of
% the filter using the vertical scroll bar. Horizontal scroll bar allows you
% to see image for each colour channel.
% All data is saved in vf_ISSPlotObject which is global object

%% Set up figure for colour channel 1
%Global paramater as plot has multiple uicontrol, makes it easier to share data between them.
global vf_ISSPlotObject

name = ['Tile ' num2str(vf_ISSPlotObject.t) ', Round ' num2str(vf_ISSPlotObject.r)];
vf_ISSPlotObject.fh = figure('units','pixels','position',...
    [500 200 800 600],'name',name,'numbertitle','off');  %Left, Bottom, Width, Height

vf_ISSPlotObject.c = 1;
vf_ISSPlotObject.Filter = 0;    %Start with unfiltered image
filter_and_plot();


%% uicontrol
%Slider to change colour channel
vf_ISSPlotObject.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',vf_ISSPlotObject.nChannels,'val',1,...
                 'sliderstep',[1/(vf_ISSPlotObject.nChannels-1) 1/(vf_ISSPlotObject.nChannels-1)],...
                 'Tag','slider1','callback',{@sl_call}); 
             
%Button to toggle filtering            
vf_ISSPlotObject.sl2 = uicontrol('style','togglebutton',...
                 'unit','pix',...
                 'position',[20 40 60 20],...
                 'val',vf_ISSPlotObject.Filter,...
                 'String','Filter',...
                 'Tag','button1','callback',{@sl2_call});             

%Vertical scroll bar to change filter radius
%Use JSlider as allows tick labels.             
jSlider = javax.swing.JSlider;         
[jhSlider, hContainer] = javacomponent(jSlider);

%As filter for Dapi channel is different, need to change slider range.
if vf_ISSPlotObject.c == vf_ISSPlotObject.DapiChannel && vf_ISSPlotObject.r == vf_ISSPlotObject.AnchorRound
        set(jSlider, 'Value',vf_ISSPlotObject.DapiR, 'Orientation',jSlider.VERTICAL, 'MajorTickSpacing',5,'PaintTicks',true,...
    'PaintLabels',true,'maximum',95,'minimum',5,'SnapToTicks',true,'MinorTickSpacing',1);
else
    set(jSlider, 'Value',vf_ISSPlotObject.R1, 'Orientation',jSlider.VERTICAL, 'MajorTickSpacing',1,'PaintTicks',true,...
    'PaintLabels',true,'maximum',20,'minimum',1,'SnapToTicks',true);
end

set(hContainer,'position',[740,50,40,400]);
hjSlider = handle(jSlider, 'CallbackProperties');
set(hjSlider, 'StateChangedCallback', @sl3_call);
vf_ISSPlotObject.jSlider = jSlider;

%Slider title
box_title = 'Filter Radius';
box = uicontrol('Style', 'text', 'Position', [740,450,50,25]);
set(box,'string', sprintf('%s', box_title));

set( findall( vf_ISSPlotObject.fh, '-property', 'Units' ), 'Units', 'Normalized' )
drawnow;
end

function [] = sl_call(h,eventdata)
% Callback for the horizontal slider.
global vf_ISSPlotObject
vf_ISSPlotObject.c = round(h.Value);

%Filter radius has different range for Dapi channel
if vf_ISSPlotObject.r == vf_ISSPlotObject.AnchorRound
    if vf_ISSPlotObject.c == vf_ISSPlotObject.DapiChannel
        set(vf_ISSPlotObject.jSlider,'maximum',95,'minimum',5,'Value',vf_ISSPlotObject.DapiR,...
            'MajorTickSpacing',5,'MinorTickSpacing',1,'LabelTable',vf_ISSPlotObject.jSlider.createStandardLabels(5));
    else
        set(vf_ISSPlotObject.jSlider, 'MajorTickSpacing',1,'maximum',20,'minimum',1,...
            'LabelTable',vf_ISSPlotObject.jSlider.createStandardLabels(1),'Value',vf_ISSPlotObject.R1);
    end
end

filter_and_plot();
end


function [] = sl2_call(h,eventdata)
% Callback for the button.
global vf_ISSPlotObject
vf_ISSPlotObject.Filter = h.Value;

filter_and_plot();  
end

function sl3_call(h,eventdata)
% Callback for vertical slider
global vf_ISSPlotObject

%Only attempt filtering once stopped scrolling, else slow
if ~h.getValueIsAdjusting
    if vf_ISSPlotObject.r == vf_ISSPlotObject.AnchorRound && vf_ISSPlotObject.c == vf_ISSPlotObject.DapiChannel
        vf_ISSPlotObject.DapiR = round(h.getValue);
    else
        vf_ISSPlotObject.R1 = round(h.getValue);
        vf_ISSPlotObject.R2 = vf_ISSPlotObject.R1*2;    %Maybe add second slider for R2??
    end
    filter_and_plot();
end
end

function [] = filter_and_plot()
%Given updated paramaters of vf_ISSPlotObject, given by callback functions,
%this plots the image.
global vf_ISSPlotObject

%Get correct image - filtered or not
if vf_ISSPlotObject.Filter == 0
    vf_ISSPlotObject.Image = vf_ISSPlotObject.RawImages(:,:,vf_ISSPlotObject.c);
else
    if vf_ISSPlotObject.r == vf_ISSPlotObject.AnchorRound &&...
            vf_ISSPlotObject.c == vf_ISSPlotObject.DapiChannel
        SE = get_filter(vf_ISSPlotObject.DapiR);
        vf_ISSPlotObject.FilterImages(:,:,vf_ISSPlotObject.c) = ...
            imtophat(vf_ISSPlotObject.RawImages(:,:,vf_ISSPlotObject.c), SE);
    else
        SE = get_filter(vf_ISSPlotObject.R1,vf_ISSPlotObject.R2);
        I_mod = single(padarray(vf_ISSPlotObject.RawImages(:,:,vf_ISSPlotObject.c),(size(SE)-1)/2,'replicate','both'));
        vf_ISSPlotObject.FilterImages(:,:,vf_ISSPlotObject.c) = convn(I_mod,SE,'valid')*vf_ISSPlotObject.Scale;
        clearvars I_mod
    end
    vf_ISSPlotObject.Image = vf_ISSPlotObject.FilterImages(:,:,vf_ISSPlotObject.c);
end

%Plot image
vf_ISSPlotObject.Background = imagesc(vf_ISSPlotObject.Image); hold on; colormap bone;
caxis([0 max(vf_ISSPlotObject.Image(:))]);

%Title
if vf_ISSPlotObject.r == vf_ISSPlotObject.AnchorRound && vf_ISSPlotObject.c==vf_ISSPlotObject.DapiChannel
    vf_ISSPlotObject.Title = title(['Tile ' num2str(vf_ISSPlotObject.t) ', Dapi Channel']);
elseif vf_ISSPlotObject.r == vf_ISSPlotObject.AnchorRound && vf_ISSPlotObject.c==vf_ISSPlotObject.AnchorChannel
    vf_ISSPlotObject.Title = title(['Tile ' num2str(vf_ISSPlotObject.t) ', Anchor Channel']);
else
    vf_ISSPlotObject.Title = title(['Tile ' num2str(vf_ISSPlotObject.t) ', Round '...
        num2str(vf_ISSPlotObject.r) ', Colour Channel ' num2str(vf_ISSPlotObject.c)]);
end
drawnow;
end

