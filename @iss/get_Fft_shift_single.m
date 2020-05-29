function [shift, score] = get_Fft_shift_single(o,t1,r1,c1,t2,r2,c2,Method,direction)
%% [shift, score] = get_Fft_shift_single(o,t1,r1,c1,t2,r2,c2,Method,direction)
%finds initial shift through Fft method from tile t1, round r1, channel c1 to t2,r2,c2.
%Only runs if other method fails.
%Method is 'Register' or 'FindSpots'
%direction is 'South' or 'East', only relevant for Register

%%
if nargin<9 || isempty(direction)
    if strcmpi(Method, 'Register')
        error('Method is Register but have not specified direction');
    else
        direction = nan;
    end
end

Graphics = o.Graphics;
o.Graphics = 2;

FileName = o.TileFiles{r1,t1};
TifObj = Tiff(FileName);
TifObj.setDirectory(o.FirstBaseChannel + c1 - 1);
%Deliberately set all negatives to 0.
ReferenceIm1 = uint16(TifObj.read())-o.TilePixelValueShift; 
FileName = o.TileFiles{r2,t2};
TifObj = Tiff(FileName);
TifObj.setDirectory(o.FirstBaseChannel + c2 - 1);
ReferenceIm2 = uint16(TifObj.read())-o.TilePixelValueShift;


if strcmpi(Method, 'Register')
    %Set CorrThresh to 0 so definitely get answer
    [shift, score] = o.ImRegFft2_Register(ReferenceIm1,ReferenceIm2, 0, o.RegMinSize, direction);
    score = o.RegAbsoluteMinScore+score;        %Shift score so not considered outlier based on score
elseif strcmpi(Method, 'FindSpots')
    [shift, score] = o.ImRegFft2_FindSpots(ReferenceIm1,ReferenceIm2, 0, o.RegMinSize);
    score = o.FindSpotsAbsoluteMinScore+score;
end    

o.Graphics = Graphics;
end