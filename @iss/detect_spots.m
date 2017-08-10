function [PeakPos, Isolated] = detect_spots(o, Image)
% [PeaksPos, Isolated] = o.DetectSpots(Image)
% 
% find positions of spots corresponding to RNA detections in a b/w image.
% The input SHOULD ALREADY HAVE BEEN TOP-HAT FILTERED (radius 3). If you
% want to smooth the image, do that first, this function won't.
% If you want to run it on generic images run DetectSpotsSingleTime
% works by finding local maxima and keeping only those where the top hat
% filter exceeds a threshold. 
%
% 
% 
% Also finds "isolated spots" for which morphological opening is less than
% a threshold. 
% 
% both thresholds are ABSOLUTE
% 
% Parameters (in iss struct):
%  DetectionThresh: threshold for top-hat filter when detecting a spot (try 0.5 * 98th Prctile)
%  DetectionRadius: radius of top-hat filter when detecting a spot (default 3)
%  IsolationThresh: threshold for morphological opening when determining if
%    a spot is isolated (try 0.05 * 98th prctile)
%  IsolationRadius: size of opening filter (default 5)
%
% For now only works on 2d but should be extendable
%
% outputs:
% PeakPos is n by 2 array of coordinates (will be y first if plotting)
% Isolation is binary array saying which are isolated
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html


%% find peaks
% first do morphological filtering
se1 = strel('disk', o.DetectionRadius);
Dilate = imdilate(Image, se1);

% local maxima are where image=dilation
Small = 1e-6; % just a small number, for computing local maxima: shouldn't matter what it is
MaxPixels = find(Image + Small >= Dilate & Image>o.DetectionThresh); 

[yPeak, xPeak] = ind2sub(size(Image), MaxPixels);
PeakPos = [yPeak, xPeak];

%% now find isolated peaks by annular filtering
if nargout==1; return; end

% first make annular filter
[xr, yr] = meshgrid(-o.IsolationRadius2:o.IsolationRadius2);
Annulus = (xr.^2 + yr.^2)<=o.IsolationRadius2.^2 & (xr.^2 + yr.^2)>o.IsolationRadius1.^2;

% filter the image with it
AnnularFiltered = imfilter(Image, double(Annulus)/sum(Annulus(:)));

% now threshold
%ScaledIsolationThresh = Range * [1-o.IsolationThresh; o.IsolationThresh];
Isolated = (AnnularFiltered(MaxPixels) < o.IsolationThresh);

% now  plot detected points on top of original image
if o.Graphics==2
    figure(50965467); clf; 
    imagesc(Image); hold on; colormap hot
    plot(xPeak(Isolated), yPeak(Isolated), 'gx');
    plot(xPeak(~Isolated), yPeak(~Isolated), 'wx');
    legend('Isolated', 'Not isolated');
    drawnow
end

return

% plot(xPeakBest, yPeakBest, 'wx');