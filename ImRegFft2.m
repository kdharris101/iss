function [shift, cc, fa1, fa2] = ImRegFft2(Im1, Im2, CorrThresh, MinSize)
% [shift, cc, f1, ft2] = ImRegFft2(Im1, Im2, CorrThresh)
%
% do image registration via fft convolution, finding match as point of 
% maximum correlation in the unwhitened images
%
% If Im1 and Im2 both come from the same global image, shift
% is the position of Im2's origin - Im1's origin.
%
% Equivalently, shift the vector such that Im2(x-shift) = Im1(x)
% and Im2(x) = Im1(x+shift) [approximately]
%
% No match if correl<CorrThresh, returns [nan nan].
%
% Correlation returned as cc. NOTE if you pass a 2-element
% vector to CorrThresh, the second entry is an extra-stringent threshold it
% uses for offsets of exactly [0 0], which is often obtained spuriously in
% microscope images. Default is [.2 .6]
%
% MinSize is a number of pixels you need to have matching before it can
% give you a good score (used to regularize the correlation)
%
% If instead of a matrix you pass a 2-element cell array for Im1 or Im2,
% this contains the fft and energy arrays, to save time. These are
% optionally returned as fa1 and fa2.
% 
% Kenneth D. Harris, 9/8/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 

% not tapering images yet but could


if nargin<3; CorrThresh = [.2 .6]; end
if length(CorrThresh)<2; CorrThresh = CorrThresh*[1, 1]; end

if nargin<4
    MinSize = 100;
end

nTries = 13; % how many local maxima to try for CorrThresh before giving up

if iscell(Im1)
    sz = size(Im1{1}, 1)/2;% /2 because it was zero-padded
else
    sz = size(Im1,1); 
end

%%
if ~iscell(Im1)
    % convert to double because matlab has all sorts of problems with integer data types
    Im1d = double(Im1);

    % create arrays of z-scored original images
    Im1z = (Im1d - mean(Im1d(:)))/std(Im1d(:));

    % zero pad them
    Im1zp = zeros(sz*2);
    Im1zp(1:sz,1:sz) = Im1z;
    
    % Fourier 
    f1 = fft2(Im1zp);

    % compute total energy in sub-images of different sizes
    % first make indefinite integrals of energy, starting with a zero:
    Cum1 = zeros(sz+1,sz+1);
    Cum1(2:sz+1,2:sz+1) = cumsum(cumsum(Im1z.^2,1),2);

    % next find box edges (inclusive), as a function of dy and dx. 0 or sz+1 means
    % no overlap
    Box1Top = [1:sz, ones(1,sz)]';
    Box1Bot = [sz*ones(1,sz) , 0:(sz-1)]';
    Box1Left = [1:sz, ones(1,sz)];
    Box1Right = [sz*ones(1,sz) , (0:sz-1)];

    % finally, doing the 2d definite integral means a difference of a difference
    Energy1 = Cum1(Box1Top,Box1Left) + Cum1(Box1Bot+1,Box1Right+1)...
            - Cum1(Box1Top,Box1Right+1) - Cum1(Box1Bot+1,Box1Left);
else
    f1 = Im1{1};
    Energy1 = Im1{2};
end


% now for image 2 - note box computation is different.
if ~iscell(Im2)
    Im2d = double(Im2);
    Im2z = (Im2d - mean(Im2d(:)))/std(Im2d(:));
    Im2zp = zeros(sz*2);
    Im2zp(1:sz,1:sz) = Im2z;
    f2 = fft2(Im2zp);

    Cum2 = zeros(sz+1,sz+1);
    Cum2(2:sz+1,2:sz+1) = cumsum(cumsum(Im2z.^2,1),2);
    Box2Top = [ones(1,sz), (sz+1):-1:2]';
    Box2Bot = [sz:-1:1, sz*ones(1,sz)]';
    Box2Left = [ones(1,sz), (sz+1):-1:2];
    Box2Right = [sz:-1:1, sz*ones(1,sz)];

    Energy2 = Cum2(Box2Bot+1,Box2Right+1) + Cum2(Box2Top,Box2Left) ...
            - Cum2(Box2Top,Box2Right+1) - Cum2(Box2Bot+1,Box2Left);
else
    f2 = Im2{1};
    Energy2 = Im2{2};
end

% convolve
Conv = ifft2(f1 .* conj(f2));

% compute correlation for each shift
Correl = (Conv./(MinSize + sqrt(Energy1.*Energy2)));

[cc, MaxShift] = max(Correl(:));

% if found zero shift, did you pass the stringent threshold?
if MaxShift==1
    if cc>=CorrThresh(2)
        [dy0, dx0] = ind2sub(size(Conv), MaxShift);
        shift = mod([dy0, dx0] +sz, sz*2) - sz - 1;  
    else
        % try second best
        [sorted, order] = sort(Correl(:), 'descend');
        cc = sorted(2);
        MaxShift = order(2);
    end
end

if MaxShift~=1  % including if you just avoided the top one
    if cc>CorrThresh(1)
        [dy0, dx0] = ind2sub(size(Conv), MaxShift);
        shift = mod([dy0, dx0] +sz, sz*2) - sz - 1;  
    else
        shift = [NaN, NaN];
    end
end

% optional pre-computation outputs:
if nargout>=3
    fa1 = {f1; Energy1};
end

if nargout>=4
    fa2 = {f2; Energy2};
end
