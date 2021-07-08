function SE = get_filter(R1,R2)
%% SE = get_filter(R1,R2)
%
% function that produces structuring element for filtering. If one
% argument, gives Dapi filter, otherwise gives imaging round / 
% anchor channel filter.
% Between r=R1 and r=R2, filter is negative. For r<R1, filter is positive.
% Overall, filter sums to 1.
%
% R1: inner radius of filter
% R2: outer radius of filter (usually 2*R1). 
% SE: filter

%%
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