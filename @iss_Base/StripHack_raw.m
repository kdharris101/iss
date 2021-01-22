function [RawIm,ChangeColumns] = StripHack_raw(o,RawIm)
%% o = o.StripHack_raw(RawIm);
%   finds all columns in RawIm where each row is identical and then sets
%   this column to the nearest normal column. Basically 'repeat padding'.
%   o: iss object
%   RawIm: the fstack_modified image, before filtering.
ChangeColumns = find(std(double(RawIm))==0);
GoodColumns = setdiff(1:o.TileSz,ChangeColumns);
for col = ChangeColumns
    [~,NearestGoodColIndex] = min(abs(GoodColumns-col));
    NearestGoodCol = GoodColumns(NearestGoodColIndex);
    RawIm(:,col) = RawIm(:,NearestGoodCol);
end
end

