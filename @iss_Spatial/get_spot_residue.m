function [OrigResidue,NewResidue] = ...
    get_spot_residue(o,LocalYX,SpotColors,LocalPeakInd,PeakCodeNo,...
    PeakCoef,BackgroundCodeNo)
%% [OrigResidue,NewResidue] = ...
% get_spot_residue(o,LocalYX,SpotColors,LocalPeakInd,PeakCodeNo,...
% PeakCoef,BackgroundCodeNo)
% Update Spot Colors by for each peak, subtracting SpotShapeDiamxSpotShapeDiam
% structuring element around it. Also update background in this area. 
% Then return residual of this area before and after. 
% o: iss object
% LocalYX: the YX location of all pixels in the tile.
% SpotColors: SpotColors(s,b,r) is the z-scored intensity in channel b,
% round r at LocalYX(s,:).
% Outputs of get_conv_peaks:
% LocalPeakInd: The peaks found are located at LocalYX(LocalPeakInd,:).
% PeakCodeNo: The peaks found correspond to o.SpatialBledCodes(PeakCodeNo,:,:).
% PeakCoef: PeakCoef(p,n) is the value of the convolution for peak p with
% spot shape eigenvector n i.e. o.SpotShapeEigenvectors(:,n).
% BackgroundCodeNo: o.SpatialBledCodes(BackgroundCodeNo,:,:) are the
% background vectors.

% Lots of confusing indexes here
nPeaks = length(PeakCodeNo);
%Get YX locations and index in image first for speed.
YXLoc = zeros(o.SpotShapeDiam^2*nPeaks,2);
for s=1:nPeaks   
    YRange = LocalYX(LocalPeakInd(s),2)-o.SpotShapeRad:LocalYX(LocalPeakInd(s),2)+o.SpotShapeRad;
    XRange = LocalYX(LocalPeakInd(s),1)-o.SpotShapeRad:LocalYX(LocalPeakInd(s),1)+o.SpotShapeRad;
    [A,B] = meshgrid(YRange,XRange);
    c=cat(2,A',B');
    YXLoc((s-1)*o.SpotShapeDiam^2+1:s*o.SpotShapeDiam^2,:) = reshape(c,[],2);      %YX pos of all NxN pixels centered on peak.  
end

YXLoc = [YXLoc(:,2),YXLoc(:,1)];    %Have to flip for indexing  
[~,LocalIndexFull]=ismember(YXLoc,LocalYX,'rows');
SpotShapes = reshape(o.SpotShapeEigenvectors(:,1:o.nShapeEigUse),...
        [o.SpotShapeDiam,o.SpotShapeDiam,o.nShapeEigUse]);
%Because conv rotates one of images so not just cross multilying:
SpotShapesRot = zeros(size(SpotShapes));
for i=1:o.nShapeEigUse
    SpotShapesRot(:,:,i) = imrotate(SpotShapes(:,:,i)',180);
end
PeakCoef = reshape(PeakCoef,nPeaks,1,o.nShapeEigUse);
OrigResidue = zeros(nPeaks,o.SpotShapeDiam,o.SpotShapeDiam);
NewResidue = zeros(nPeaks,o.SpotShapeDiam,o.SpotShapeDiam);
for s=1:nPeaks
    LocalIndex = LocalIndexFull((s-1)*o.SpotShapeDiam^2+1:s*o.SpotShapeDiam^2);
    PredCode = sum(PeakCoef(s,:,:).*SpotShapesRot,3);
    PredCode(LocalIndex==0)=0;
    OrigCode = zeros(o.SpotShapeDiam^2,o.nBP,o.nRounds);
    OrigCode(LocalIndex~=0,:,:) = SpotColors(LocalIndex(LocalIndex~=0),:,:);
    OrigResidue(s,:,:) = reshape(vecnorm(OrigCode(:,:),2,2),...
        [o.SpotShapeDiam,o.SpotShapeDiam]);
    NewCode = zeros(o.SpotShapeDiam^2,o.nBP,o.nRounds);
    NewCode(LocalIndex~=0,:,:) = SpotColors(LocalIndex(LocalIndex~=0),:,:)-...
        PredCode(LocalIndex~=0).*o.SpatialBledCodes(PeakCodeNo(s),:,:);  
    NewCode2 = update_spot_background(o,NewCode,BackgroundCodeNo);
    NewCode(LocalIndex~=0,:,:) = NewCode2(LocalIndex~=0,:,:);
    NewResidue(s,:,:) = reshape(vecnorm(NewCode(:,:),2,2),...
        [o.SpotShapeDiam,o.SpotShapeDiam]);
end

end
