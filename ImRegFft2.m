function [shift, cc] = ImRegFft2(Im1, Im2, CorrThresh, MinSize)
% shift = ImRegFft2(Im1, Im2, CorrThresh)
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
% Kenneth D. Harris, 9/8/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
 

% not tapering images yet but could


if nargin<3; CorrThresh = [.2 .6]; end
if length(CorrThresh)<2; CorrThresh = CorrThresh*[1, 1]; end

if nargin<4
    MinSize = 100;
end

nTries = 13; % how many local maxima to try for CorrThresh before giving up

sz = size(Im1, 1);

%%
% convert to double because matlab has all sorts of problems with integer data types
Im1d = double(Im1);
Im2d = double(Im2);

% create arrays of z-scored original images
Im1z = (Im1d - mean(Im1d(:)))/std(Im1d(:));
Im2z = (Im2d - mean(Im2d(:)))/std(Im2d(:));

% zero pad them
Im1zp = zeros(sz*2);
Im2zp = zeros(sz*2);
Im1zp(1:sz,1:sz) = Im1z;
Im2zp(1:sz,1:sz) = Im2z;

% compute total energy in sub-images of different sizes
% first make indefinite integrals of energy, starting with a zero:
Cum1 = zeros(sz+1,sz+1);
Cum2 = zeros(sz+1,sz+1);
Cum1(2:sz+1,2:sz+1) = cumsum(cumsum(Im1z.^2,1),2);
Cum2(2:sz+1,2:sz+1) = cumsum(cumsum(Im2z.^2,1),2);

if 1 % THIS WAY IS FASTER!
    % now find box edges (inclusive), as a function of dy and dx. 0 or sz+1 means
    % no overlap
    Box1Top = [1:sz, ones(1,sz)]';
    Box1Bot = [sz*ones(1,sz) , 0:(sz-1)]';
    Box1Left = [1:sz, ones(1,sz)];
    Box1Right = [sz*ones(1,sz) , (0:sz-1)];

    % doing the 2d definite integral means a difference of a difference
        
    Energy1 = Cum1(Box1Top,Box1Left) + Cum1(Box1Bot+1,Box1Right+1)...
            - Cum1(Box1Top,Box1Right+1) - Cum1(Box1Bot+1,Box1Left);

    Box2Top = [ones(1,sz), (sz+1):-1:2]';
    Box2Bot = [sz:-1:1, sz*ones(1,sz)]';
    Box2Left = [ones(1,sz), (sz+1):-1:2];
    Box2Right = [sz:-1:1, sz*ones(1,sz)];

    Energy2 = Cum2(Box2Bot+1,Box2Right+1) + Cum2(Box2Top,Box2Left) ...
            - Cum2(Box2Top,Box2Right+1) - Cum2(Box2Bot+1,Box2Left);
else % THIS WAY IS SLOWER AND STUPIDER
    % find energy in boxes. need a different algorithm for each quadrant
    % sum of energy over square [a,b]->[c,d] inclusive is 
    % Cum1(a,b) + Cum1(c+1,d+1) -Cum1(a,d+1) - Cum1(c+1,b)
    
    % box1
    
    % top left quad: Box1 goes from (y,x) to (sz,sz)
    
    Energy1(1:sz, 1:sz) = Cum1(1:sz,1:sz) + Cum1(sz+1,sz+1) ...
                        - Cum1(1:sz, sz+1) - Cum1(sz+1, 1:sz); 
    % Box2 from (1,1) to (sz-y+1, sz-x+1)
    Energy2(1:sz, 1:sz) = Cum2(1,1) + Cum2(sz+1:-1:2,sz+1:-1:2) ...
                        - Cum2(1,sz+1:-1:2) - Cum2(sz+1:-1:2, 1);
    
    % top right quad: Box1 goes from (y,1) to (sz,x-1)
    Energy1(1:sz, sz+1:2*sz) = Cum1(1:sz,1) + Cum1(sz+1,1:sz) ...
                             - Cum1(1:sz,1:sz) - Cum1(sz+1,1); 
    % Box2 from (1,sz+2-x) to (sz-y+1,sz)
    Energy2(1:sz, sz+1:2*sz) = Cum2(1,sz+1:-1:2) + Cum2(sz+1:-1:2,sz+1) ...
                             - Cum2(1,sz+1) - Cum2(sz+1:-1:2, sz+1:-1:2);
                         
	% bottom left quad: Box1 goes from (1,x) to (y-1,sz)
    Energy1(sz+1:2*sz, 1:sz) = Cum1(1,1:sz) + Cum1(1:sz,sz+1) ...
                             - Cum1(1:sz, 1:sz) - Cum1(1, sz+1); 
    % Box2 from (sz+2-y,1) to (sz, sz-x+1)
    Energy2(sz+1:2*sz, 1:sz) = Cum2(sz+1:-1:2,1) + Cum2(sz+1,sz+1:-1:2) ...
                             - Cum2(sz+1,1) - Cum2(sz+1:-1:2, sz+1:-1:2);

    % bottom right quad: Box1 goes from (1,1) to (y-1,x-1)
    Energy1(sz+1:2*sz, sz+1:2*sz) = Cum1(1,1) + Cum1(1:sz,1:sz) ...
                             - Cum1(1,1:sz) - Cum1(1:sz, 1); 
    % Box2 from (sz+2-y,sz+2-x) to (sz, sz)
    Energy2(sz+1:2*sz, sz+1:2*sz) = Cum2(sz+1:-1:2,sz+1:-1:2) + Cum2(sz+1,sz+1) ...
                                  - Cum2(sz+1:-1:2,sz+1) - Cum2(sz+1, sz+1:-1:2);

end

% Fourier and convolve
f1 = fft2(Im1zp);
f2 = fft2(Im2zp);
Conv = ifft2(f1 .* conj(f2));

% compute correlation for each shift
Correl = (Conv./(MinSize + sqrt(Energy1.*Energy2)));

[cc, MaxShift] = max(Correl(:));

% set shift NaN if below threshold
if cc < CorrThresh(1) || (MaxShift==1 && cc<CorrThresh(2))
    shift = [NaN, NaN];
else
    [dy0, dx0] = ind2sub(size(Conv), MaxShift);
    shift = mod([dy0, dx0] +sz, sz*2) - sz - 1;    
end
%%
return


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