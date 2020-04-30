function [shift, score] = get_Fft_shift_single(o,t1,t2,direction)
%finds initial shift through Fft method. Only runs if other method fails.

Graphics = o.Graphics;
o.Graphics = 2;

FileName = o.TileFiles{o.ReferenceRound,t1};
TifObj = Tiff(FileName);
TifObj.setDirectory(o.FirstBaseChannel + o.ReferenceChannel - 1);
%Deliberately set all negatives to 0.
ReferenceIm1 = uint16(TifObj.read())-o.TilePixelValueShift; 
FileName = o.TileFiles{o.ReferenceRound,t2};
TifObj = Tiff(FileName);
TifObj.setDirectory(o.FirstBaseChannel + o.ReferenceChannel - 1);
ReferenceIm2 = uint16(TifObj.read())-o.TilePixelValueShift;

%Set CorrThresh to 0 so definitely get answer
[shift, score] = o.ImRegFft2_Register(ReferenceIm1,...
            ReferenceIm2, o.RegCorrThresh*0, o.RegMinSize, direction);
score = o.RegAbsoluteMinScore+score;        %Shift score so not considered outlier based on score
o.Graphics = Graphics;
end