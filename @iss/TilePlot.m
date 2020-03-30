function [] = TilePlot(o,FilterImage,RawImage,r,nChannels)
% Used for debugging in extract_and_filter. Plots result of filtering and
% auto scaling for tile 1 of round 1 or anchor round.
% You can only use the buttons and scrollers when whole thing finished
% though.
name = ['Tile 1, Round ' num2str(r)];
S.fh = figure('units','pixels','position',[500 200 800 600],'name',name,'numbertitle','off');  %Left, Bottom, Width, Height
S.FilterImage = FilterImage;
S.RawImage = RawImage;

S.r = r;
S.ReferenceRound = o.ReferenceRound;
S.DapiChannel = o.DapiChannel;
S.AnchorChannel = o.AnchorChannel;
S.Channel = 1;
S.Image = S.FilterImage;
S.Background = imagesc(S.Image(:,:,S.Channel)); hold on; colormap bone;
caxis([0 max(max(S.Image(:,:,S.Channel)))]);
if S.r == S.ReferenceRound && S.DapiChannel==1
    S.Title = title('Tile 1, Dapi Channel');
elseif S.r == S.ReferenceRound && S.AnchorChannel==1
    S.Title = title('Tile 1, Anchor Channel');
else
    S.Title = title(['Tile 1, Round ' num2str(S.r) ', Colour Channel 1']);
end
xlim([0 o.TileSz]);
ylim([0 o.TileSz]);

S.SliderName = ['slider',num2str(r)];   
S.ButtonName = ['button',num2str(r)]; 
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',nChannels,'val',1,...
                 'sliderstep',[1/(nChannels-1) 1/(nChannels-1)],...
                 'Tag',S.SliderName,'callback',{@sl_call,S}); 
             
            
S.sl2 = uicontrol('style','togglebutton',...
                 'unit','pix',...
                 'position',[20 40 60 20],...
                 'val',1,...
                 'String','Filter',...
                 'Tag',S.ButtonName,'callback',{@sl2_call,S});             
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )
drawnow;
end

function [] = sl_call(h,eventdata,S)
% Callback for the slider.
%[h,S] = varargin{[1,3]};  % calling handle and data structure.
Value = h.Value;
h2 = findobj('Tag',S.ButtonName);
if ~isempty(h2.UserData)
    h = h2;
    S = h.UserData;
else
    h2 = findobj('Tag',S.SliderName);
    if ~isempty(h2.UserData)
        h = h2;
        S = h.UserData;
    end
end

S.Channel = round(Value);
S.Background = imagesc(S.Image(:,:,S.Channel)); hold on; colormap bone;
caxis([0 max(max(S.Image(:,:,S.Channel)))]);
if S.r == S.ReferenceRound && S.DapiChannel==S.Channel
    S.Title = title('Tile 1, Dapi Channel');
elseif S.r == S.ReferenceRound && S.AnchorChannel==S.Channel
    S.Title = title('Tile 1, Anchor Channel');
else
    S.Title = title(['Tile 1, Round ' num2str(S.r) ', Colour Channel ' num2str(S.Channel)]);
end
h.UserData = S;
drawnow;
end

function [] = sl2_call(h,eventdata,S)
% Callback for the slider.
%[h,S] = varargin{[1,3]};  % calling handle and data structure.
Value = h.Value;
h2 = findobj('Tag',S.SliderName);
if ~isempty(h2.UserData)
    h = h2;
    S = h.UserData;
else
    h2 = findobj('Tag',S.ButtonName);
    if ~isempty(h2.UserData)
        h = h2;
        S = h.UserData;
    end
end
    
if Value==1
    S.Image = S.FilterImage;
else
    S.Image = S.RawImage;
end
S.Background = imagesc(S.Image(:,:,S.Channel)); hold on; colormap bone;
caxis([0 max(max(S.Image(:,:,S.Channel)))]);
if S.r == S.ReferenceRound && S.DapiChannel==S.Channel
    S.Title = title('Tile 1, Dapi Channel');
elseif S.r == S.ReferenceRound && S.AnchorChannel==S.Channel
    S.Title = title('Tile 1, Anchor Channel');
else
    S.Title = title(['Tile 1, Round ' num2str(S.r) ', Colour Channel ' num2str(S.Channel)]);
end
h.UserData = S;
drawnow;
end

