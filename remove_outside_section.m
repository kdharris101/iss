function [ImMod,PixelLabel] = remove_outside_section(Im,GradSz,SmoothSz)
% Idea is to remove part of image outside section by noting that this part
% of image has a much lower range of intensity values than the part within
% the section.
% Im: Image for which, we want to extract section.
% RangeSz: Filter size of rangefilt.
% CloseSz: Filter size of closefilt.
% ImMod: Im where pixels such that PixelLabel=true are set to 0.

if nargin<2 || isempty(GradSz)
    GradSz = 5;
end
if nargin<3 || isempty(SmoothSz)
    SmoothSz = 200;
end
Im = double(Im);
AvSz = 100;
AvSE = ones(AvSz)/(AvSz*AvSz);
GradSE = repelem([-1,0,1],GradSz);
GradSE = repmat(GradSE,GradSz,1);
GradSE = GradSE/sum(abs(GradSE(:)));
ImAverage = imfilter(Im,AvSE,'same','replicate');
ImGrad = imfilter(ImAverage,GradSE,'same','replicate');
ImAvThresh = repelem(prctile(ImAverage',90)',1,size(Im,2));
ImGradThresh = repelem(prctile(ImGrad',90)',1,size(Im,2));
Use = Im>ImAvThresh | ImGrad>ImGradThresh;
if sum(Use(:))/length(Use(:))<0.05
    %If less than 5% of pixels to be changed then don't change any
    PixelLabel = false(size(Im));
    ImMod = Im;
else
    [a,b] = ind2sub(size(Im),find(Use));
    PixelLabel = false(size(Im));
    for i=1:size(Im,1)
        PixelLabel(i,1:min(b(a==i))) = true;
    end
    
    %Smooth PixelLabel    
    SmoothSE = ones(SmoothSz);
    SmoothSE = SmoothSE/(SmoothSz*SmoothSz);
    PixelLabel = imfilter(PixelLabel,SmoothSE,'same','replicate');
    ImMod = Im;
    ImMod(PixelLabel) = 0;
end
% figure; imagesc(PixelLabel);
% figure; imagesc(ImMod);


% figure; imagesc(ImMod);
% hold on
% scatter(Seed_x,Seed_y,100,'rx');
% hold off
% set(gca, 'YDir', 'normal');


end

