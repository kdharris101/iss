function SpotColors = ...
    update_spot_colors(o,LocalYX,SpotColors,LocalPeakInd,PeakCodeNo,PeakCoef)
%% SpotColors = ...
% o.update_spot_colors(SpotColors,LocalYX,LocalPeakInd,PeakCodeNo,PeakCoef)
% Update Spot Colors by for each peak, subtracting SpotShapeDiamxSpotShapeDiam
% structuring element around it.
% o: iss object
% LocalYX: the YX location of all pixels in the tile.
% SpotColors: SpotColors(s,b,r) is the z-scored intensity in channel b,
% round r at LocalYX(s,:).
% Outputs of get_conv_peaks:
% LocalPeakInd: The peaks found are located at LocalYX(LocalPeakInd,:).
% PeakCodeNo: The peaks found correspond to o.SpatialBledCodes(PeakCodeNo,:,:).
% PeakCoef: PeakCoef(p,n) is the value of the convolution for peak p with
% spot shape eigenvector n i.e. o.spShapeEigenvectors(:,n).
%% 
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
for s=1:nPeaks
    LocalIndex = LocalPeakFullShapeIndices((s-1)*o.spShapeDiam^2+1:s*o.spShapeDiam^2);
    PredCode = sum(PeakCoef(s,:,:).*SpotShapesRot,3);
    PredCode(LocalIndex==0)=0;
    NewCode = zeros(o.spShapeDiam^2,o.nBP,o.nRounds);
    NewCode(LocalIndex~=0,:,:) = SpotColors(LocalIndex(LocalIndex~=0),:,:)-...
        PredCode(LocalIndex~=0).*spBledCodes(PeakCodeNo(s),:,:);
    SpotColors(LocalIndex(LocalIndex~=0),:,:) = NewCode(LocalIndex~=0,:,:);
    % For help debugging:
    %A = SpotColors(LocalIndex,:)*o.SpatialBledCodes(PeakCodeNo(s),:)';
    %A_coefs = convn(reshape(A,[17,17]),permute(SpotShapes,[2,1,3]));  
    %A_coefs = squeeze(A_coefs(o.spShapeDiam,o.spShapeDiam,:));  %THIS SHOULD BE PeakCoef(s,:)
    %A2 = NewCode(:,:)*o.SpatialBledCodes(PeakCodeNo(s),:)';
    %A2_coefs = convn(reshape(A2,[17,17]),permute(SpotShapes,[2,1,3]));
    %A2_coefs = squeeze(A2_coefs(o.spShapeDiam,o.spShapeDiam,:)); %THIS SHOULD BE ZERO
    %Old debugging:
    %PeakSpotColor = zeros(o.spShapeDiam^2,7,7);
    %PeakSpotColor(LocalIndex~=0,:,:) = SpotColors(LocalIndex(LocalIndex~=0),:,:);
    % PeakCoefs = zeros(o.spShapeDiam^2,1);
    % for s2=1:o.spShapeDiam^2
    %     if GoodAnchorIndex(s2)~=0
    %         PeakCoefs(s2,:) = omp_free_background(o.SpatialBledCodes(BackgroundType,:)',PeakSpotColor(s2,:)',...
    %             1)';
    %     end
    % end
    % PeakCoefs = reshape(PeakCoefs,[o.spShapeDiam,o.spShapeDiam])';
   
end

end
