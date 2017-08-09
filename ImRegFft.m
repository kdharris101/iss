function [shift, cc] = ImRegFft(Im1, Im2, Direction, CorrThresh, MinSize)
% shift = ImRegFft(Im1, Im2, Direction, CorrThresh)
%
% do image registration via fft convolution. 
% assumes square inputs
%
% If Im1 and Im2 both come from the same global image, shift
% is the position of Im2's origin - Im1's origin.
%
% Equivalently, shift the vector such that Im2(x-shift) = Im1(x)
% and Im2(x) = Im1(x+shift) [approximately]
%
% If no match, returns [nan nan].
%
%
% Direction is 'se', 'nw', 'n', etc. or '' to indicate where Im2 is going to be
% relative to Im1. (Necessary because ffts work circularly.)
% n means smaller values of y; e means larger values of x. '' means close
% to centered
%
% CorrThresh (default .6): correlation of images on overlap region needs to be this high
% or no match. Correlation returned as cc. NOTE if you pass a 2-element
% vector to CorrThresh, the second entry is an extra-stringent threshold it
% uses for offsets of exactly [0 0], which is often obtained spuriously in
% microscope images.
% 
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 

% not tapering images yet but could


Interactive = 0;

if nargin<4
    CorrThresh = .6;
end

if nargin<5
    MinSize = 100;
end

nTries = 13; % how many local maxima to try for CorrThresh before giving up


% because matlab has all sorts of problems with integer data types
Im1 = double(Im1);
Im2 = double(Im2); 

sz = size(Im1, 1);

f1 = fft2(Im1);
f2 = fft2(Im2);

% to do phase correlation, whiten signals:
WhitenReg = 1e0;
wf1 = f1 ./ (abs(f1) + WhitenReg*std(Im1(:)));
wf2 = f2 ./ (abs(f2) + WhitenReg*std(Im2(:)));

Conv = ifft2(wf1 .* conj(wf2));

[~, order] = sort(Conv(:), 'descend');

shift = [nan nan]; % what to return if you find no good fits
cc = 0;
for i=1:nTries
    [dy0, dx0] = ind2sub(size(Conv), order(i));
    % have to subtract ones because matlab
    if ismember('n', Direction)
        dy = dy0-sz-1;
    elseif ismember('s', Direction)
        dy = dy0-1;
    else
        dy = mod(dy0+sz/2, sz) - sz/2 - 1;
    end
    
    if ismember('w', Direction)
        dx = dx0-sz-1;
    elseif ismember('e', Direction)
        dx = dx0-1;
    else
        dx = mod(dx0+sz/2, sz) - sz/2 - 1;
    end
    
    Im2Left = max(1,1-dx);
    Im2Right = min(sz, sz-dx);
    Im2Top = max(1,1-dy);
    Im2Bot = min(sz, sz-dy);
    
    SubIm2 = Im2(Im2Top:Im2Bot, Im2Left:Im2Right);
    SubIm1 = Im1(Im2Top+dy:Im2Bot+dy, Im2Left+dx:Im2Right+dx);

    if numel(SubIm1)==0 || numel(SubIm2)==0
        continue;
    end
    
    cc = corr(SubIm1(:), SubIm2(:));

    if Interactive
         figure(21789)
         subplot(3,1,1); imagesc(Im1); colorbar;
         subplot(3,1,2); imagesc(Im2); colorbar;
         Sc = prctile([Im1(:); Im2(:)], 100);
         subplot(3,1,3); image(cat(3, SubIm1/Sc, SubIm2/Sc, SubIm1*0));

        fprintf('shift %d %d, cc %f\n', dy, dx, cc);
        pause
    end

    
    if min(size(SubIm1))<MinSize || min(size(SubIm2))<MinSize
        continue;
    end

    if dx~=0 || dy~=0 
        % lenient first threshold for shifts not zero
        if cc>CorrThresh(1)
            shift = [dy dx];
            break;
        end
    else 
        % tougher threshold for shifts of zero
        if cc>max(CorrThresh(:))
            shift = [dy dx];
            break;
        end
    end
end
return