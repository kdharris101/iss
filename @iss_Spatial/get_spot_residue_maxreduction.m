function MaxResidueReduction = ...
    get_spot_residue_maxreduction(o,LocalYX,SpotColors,LocalPeakInd,PeakCodeNo,...
    PeakCoef,BackgroundCodeNo)
%% MaxResidueReduction = ...
% get_spot_residue(o,LocalYX,SpotColors,LocalPeakInd,PeakCodeNo,...
% PeakCoef,BackgroundCodeNo)
% Update Spot Colors by for each peak, subtracting SpotShapeDiamxSpotShapeDiam
% structuring element around it. Also update background in this area. 
% Then return value of maximum residual reduction in the area.
% o: iss object
% LocalYX: the YX location of all pixels in the tile.
% SpotColors: SpotColors(s,b,r) is the z-scored intensity in channel b,
% round r at LocalYX(s,:).
% Outputs of get_conv_peaks:
% LocalPeakInd: The peaks found are located at LocalYX(LocalPeakInd,:).
% PeakCodeNo: The peaks found correspond to o.SpatialBledCodes(PeakCodeNo,:,:).
% PeakCoef: PeakCoef(p,n) is the value of the convolution for peak p with
% spot shape eigenvector n i.e. o.spShapeEigenvectors(:,n).
% BackgroundCodeNo: o.spBledCodes(BackgroundCodeNo,:,:) are the
% background vectors.

spBledCodes = reshape(o.spBledCodes,size(o.spBledCodes,1),o.nBP,o.nRounds);
% Get LocalYX indices of all pixels in the shape around the peaks 
LocalPeakFullShapeIndices = o.get_spot_shape_indices(LocalYX,LocalPeakInd);

SpotShapes = reshape(o.spShapeEigenvectors(:,1:o.nShapeEigUse),...
        [o.spShapeDiam,o.spShapeDiam,o.nShapeEigUse]);
%Because conv rotates one of images so not just cross multilying:
SpotShapesRot = zeros(size(SpotShapes));
for i=1:o.nShapeEigUse
    SpotShapesRot(:,:,i) = imrotate(SpotShapes(:,:,i)',180);
end

nPeaks = length(LocalPeakInd);
PeakCoef = reshape(PeakCoef,nPeaks,1,o.nShapeEigUse);
OrigResidue = zeros(nPeaks,o.spShapeDiam*o.spShapeDiam);
NewResidue = zeros(nPeaks,o.spShapeDiam*o.spShapeDiam);
for s=1:nPeaks
    LocalIndex = LocalPeakFullShapeIndices((s-1)*o.spShapeDiam^2+1:s*o.spShapeDiam^2);
    PredCode = sum(PeakCoef(s,:,:).*SpotShapesRot,3);
    PredCode(LocalIndex==0)=0;
    OrigCode = zeros(o.spShapeDiam^2,o.nBP,o.nRounds);
    OrigCode(LocalIndex~=0,:,:) = SpotColors(LocalIndex(LocalIndex~=0),:,:);
    OrigResidue(s,:) = vecnorm(OrigCode(:,:),2,2);
    NewCode = zeros(o.spShapeDiam^2,o.nBP,o.nRounds);
    NewCode(LocalIndex~=0,:,:) = SpotColors(LocalIndex(LocalIndex~=0),:,:)-...
        PredCode(LocalIndex~=0).*spBledCodes(PeakCodeNo(s),:,:);  
    NewCode2 = update_spot_background(o,NewCode,BackgroundCodeNo);
    NewCode(LocalIndex~=0,:,:) = NewCode2(LocalIndex~=0,:,:);
    NewResidue(s,:) = vecnorm(NewCode(:,:),2,2);
end
MaxResidueReduction = max(OrigResidue-NewResidue,[],2);

end
