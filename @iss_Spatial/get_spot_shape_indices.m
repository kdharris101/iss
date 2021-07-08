function LocalPeakFullShapeIndices = get_spot_shape_indices(o,LocalYX,LocalPeakInd)
%% LocalPeakFullShapeIndices = o.get_spot_shape_indices(LocalYX,LocalPeakInd);
% o: iss object
% LocalYX: the YX location of all pixels in the tile.
% LocalPeakInd: The peaks found are located at LocalYX(LocalPeakInd,:).
% LocalPeakFullShapeIndices((s-1)*o.spShapeDiam^2+1:s*o.spShapeDiam^2) are
% the indices in LocalYX for all pixels in spot shape around peak at
% LocalYX(LocalPeakInd(s),:).
% LocalPeakFullShapeIndices=0 if goes off edge of tile. 
%% Lots of confusing indexes here
nPeaks = length(LocalPeakInd);

%For each peak, find xy location of all pixels surrounding it.
XYLoc = zeros(o.spShapeDiam^2*nPeaks,2);
for s=1:nPeaks   
    %For each peak, find xy location of all pixels surrounding it.
    XRange = LocalYX(LocalPeakInd(s),2)-o.spShapeRad:LocalYX(LocalPeakInd(s),2)+o.spShapeRad;
    YRange = LocalYX(LocalPeakInd(s),1)-o.spShapeRad:LocalYX(LocalPeakInd(s),1)+o.spShapeRad;
    [A,B] = meshgrid(XRange,YRange);
    c=cat(2,A',B');
    XYLoc((s-1)*o.spShapeDiam^2+1:s*o.spShapeDiam^2,:) = reshape(c,[],2);      %YX pos of all NxN pixels centered on peak.  
end
%Find index of these locations in the full tile array of YX coordinates.
YXLoc = [XYLoc(:,2),XYLoc(:,1)];
[~,LocalPeakFullShapeIndices]=ismember(YXLoc,LocalYX,'rows');
end

