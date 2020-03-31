function SE = get_filter(R1,R2)
% function that produces structuring element for filtering. If one
% argument, gives Dapi filter, otherwise gives imaging round / 
% anchor channel filter.

if nargin<2 || isempty(R2)
    SE = strel('disk', R1);
else
    h = -hanning(R2*2+1);
    h = -h/sum(h);
    h(R2+1-R1:R2+1+R1) = ...
        h(R2+1-R1:R2+1+R1)+hanning(R1*2+1)/sum(hanning(R1*2+1));
    SE = ftrans2(h');
    SE = single(SE);
end