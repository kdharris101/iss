function [shift, score] = get_Fft_shift_single(o,t1,r1,c1,t2,r2,c2,Method,direction)
%% [shift, score] = get_Fft_shift_single(o,t1,r1,c1,t2,r2,c2,Method,direction)
%finds initial shift through Fft method from tile t1, round r1, channel c1 to t2,r2,c2.
%Only runs if other method fails.
%Method is 'Register' or 'FindSpots' or 'SubPixel'
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
%Deliberately set all negatives to 0. DONT KNOW WHY I DID THIS - GOT RID OF IT
ReferenceIm1 = int32(TifObj.read())-o.TilePixelValueShift; 
FileName = o.TileFiles{r2,t2};
TifObj = Tiff(FileName);
TifObj.setDirectory(o.FirstBaseChannel + c2 - 1);
ReferenceIm2 = int32(TifObj.read())-o.TilePixelValueShift;


if strcmpi(Method, 'Register')
    %Set CorrThresh to 0 so definitely get answer
    [shift, score] = o.ImRegFft2_Register(ReferenceIm1,ReferenceIm2, 0, o.RegMinSize, direction);
    score = o.RegAbsoluteMinScore+score;        %Shift score so not considered outlier based on score
elseif strcmpi(Method, 'FindSpots')
    [shift, score] = o.ImRegFft2_FindSpots(ReferenceIm1,ReferenceIm2, 0, o.RegMinSize);
    score = o.FindSpotsAbsoluteMinScore+score;
elseif strcmpi(Method, 'SubPixel')
    %Before finding shift, need to correct for chromatic aberration by
    %scaling RefIm1 which is not in reference channel to be in anchor channel.
    %Then shift is then found as if both are in same (reference) colour channel.
    C = o.TileCentre(1);
    b=0;
    d=0;
    A1=1/o.A(c1,1);
    A2=1/o.A(c1,2);
    tform = affine2d([ ...
        A1 b 0;...
        d A2 0; ...
        C*(1-A1-b) C*(1-A2-d) 1]);
    sameAsInput = affineOutputView(size(ReferenceIm1),tform,'BoundsStyle','SameAsInput');
    ScaledReferenceIm1 = imwarp(ReferenceIm1,tform,'OutputView',sameAsInput);
    output = dftregistration(fft2(ReferenceIm1),fft2(ReferenceIm2),1000);
    score = output(1);
    shift = output(3:4);
end    

o.Graphics = Graphics;
end