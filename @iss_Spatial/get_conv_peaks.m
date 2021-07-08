function [LocalPeakInd,PeakCodeNo,PeakCoef,PeakIter,MaxResidueReduction,ConvImagesAllCodes] =...
    get_conv_peaks(o,LocalYX,SpotColors,ResidualThresh,...
    PrevLocalInd,PrevCodeNo,PrevIter,iter,BackgroundCodeNo)
%% [PeakCodeNo,PeakCoef,LocalPeakInd] =...
%    get_conv_peaks(o,LocalYX,SpotColors,CodeNumbers)
% This takes all spot colors on a tile, gets the dot product of these with
% each gene bled code. It then takes the resultant dot
% product images and convolves them with the spot shape eigenvectors.
% The location of the spots is determined by the peak in the first shape
% eigenvector convolution but coeffficients with all o.nShapeEigUse are
% returned. 
% o: iss object
% LocalYX: the YX location of all pixels in the tile.
% SpotColors: SpotColors(s,b,r) is the z-scored intensity in channel b,
% round r at LocalYX(s,:).
% ResidualThresh(s): for spot s, residual must be reduced by more than this
% for gene to be accepted. 
% PrevLocalInd: Index of peaks found on previous iterations.
% PrevCodeNo: Code number of peaks found on previous iterations.
% PrevIter: iteration where peaks found previously.
% iter: current iteration.
% LocalPeakInd: The peaks found are located at LocalYX(LocalPeakInd,:).
% PeakCodeNo: The peaks found correspond to o.SpatialBledCodes(PeakCodeNo,:,:).
% PeakCoef: PeakCoef(p,n) is the value of the convolution for peak p with
% spot shape eigenvector n i.e. o.SpotShapeEigenvectors(:,n).
% PeakIter(p) is the iteration that peak p was found on.
% MaxResidueReduction(s) is the biggest residue reduction in shape around
% spot s as a result of removing gene PeakCodeNo(s). This only includes
% spots found on this iteration. 
% ConvImagesAllCodes(:,:,g) is the convolution image for gene g with the
% first spatial eigenvector.


%% Get convolution of each code image with each shape eigenvector
%Only consider convolutions with each gene, background is considered on by
%pixel basis. 
nCodes = length(o.CharCodes);   
CodeNumbers = 1:nCodes;

SpotShapes = reshape(o.spShapeEigenvectors(:,1:o.nShapeEigUse),...
        [o.spShapeDiam,o.spShapeDiam,o.nShapeEigUse]);
ConvImages = zeros([max(LocalYX),o.nShapeEigUse,nCodes]);
for CodeNo=1:nCodes
    CodeDotProduct = SpotColors(:,:)*o.spBledCodes(CodeNumbers(CodeNo),:)';   %SLOWEST LINE - TRY TO QUICKEN BY DOING FOR ALL GENES AT ONCE 
    %Turn coefficients into an image
    CodeIm = zeros(max(LocalYX));
    Ind = sub2ind(size(CodeIm),LocalYX(:,1),LocalYX(:,2));
    CodeIm(Ind) = CodeDotProduct;
    %Convolve with each eigenvector
    ConvIm = convn(CodeIm,SpotShapes);
    %Padding doesn't work before conv so have to do after:
    ConvImages(:,:,:,CodeNo) = ConvIm(o.spShapeRad+1:end-o.spShapeRad,...
        o.spShapeRad+1:end-o.spShapeRad,:);
end

clear CodeDotProduct
%Only consider first shape as this is centered on centre pixel.
%If consider all to find peak, then spot locations weird. 
AbsConvImages = abs(squeeze(ConvImages(:,:,1,:)));

if nargout>=5
    %Find sign of first shape eigenvector
    [~,SignIndex] = max(abs(o.spShapeEigenvectors(:,1)));
    ConvSign = sign(o.spShapeEigenvectors(SignIndex,1));
    %For visualisation, want to multiply by this sign. 
    ConvImagesAllCodes = ConvSign*squeeze(ConvImages(:,:,1,:));
end
%Take maximum projection image, then find all peaks over threshold
[MaxProject,MaxProjectCodeNo] = max(AbsConvImages,[],3);
clear AbsConvImages
MaxProjectCodeNo = CodeNumbers(MaxProjectCodeNo);
Small = 1e-6;
%Want min separation of peaks to equal structuring element size
SE = strel('disk', o.spShapeRad);     
Dilate = imdilate(MaxProject, SE);
MaxProjectPeakInd = find(MaxProject + Small >= Dilate);
[~,LocalPeakInd] = ismember(MaxProjectPeakInd,Ind);
%Only keep peaks that we have a spot color for.
MaxProjectPeakInd = MaxProjectPeakInd(LocalPeakInd~=0);
LocalPeakInd = LocalPeakInd(LocalPeakInd~=0);

%Add previously found spots to those found on this iteration
LocalPeakInd = [LocalPeakInd;PrevLocalInd];
PeakCodeNo = MaxProjectCodeNo(MaxProjectPeakInd);
nPeaks_Iter = length(PeakCodeNo);
PeakCodeNo = [PeakCodeNo;PrevCodeNo];
nPeaks = length(PeakCodeNo);

%find the convolution value (coefficient) for each eigenvector for each spot
CoefY = repelem(LocalYX(LocalPeakInd,1),o.nShapeEigUse,1);
CoefX = repelem(LocalYX(LocalPeakInd,2),o.nShapeEigUse,1);
CoefEigNo = repmat(1:o.nShapeEigUse,1,nPeaks)';
CoefCodeNo = repelem(PeakCodeNo-min(CodeNumbers)+1,o.nShapeEigUse,1);
if isempty(CoefCodeNo)
    CoefCodeNo = CoefEigNo;     %When CoefCodeNo was empty, was wrong size.
end
ConvCoefsInd = sub2ind(size(ConvImages),CoefY,CoefX,CoefEigNo,CoefCodeNo);
ConvCoefs = ConvImages(ConvCoefsInd);
PeakCoef = reshape(ConvCoefs,o.nShapeEigUse,nPeaks)';

%Only accept new peaks with MaxResiudeReduction>ResidualThresh
MaxResidueReduction = get_spot_residue_maxreduction(o,LocalYX,SpotColors,...
    LocalPeakInd(1:nPeaks_Iter),PeakCodeNo(1:nPeaks_Iter),...
    PeakCoef(1:nPeaks_Iter,:),BackgroundCodeNo);
Keep = true(size(PeakCodeNo));
Keep(1:nPeaks_Iter) = MaxResidueReduction>ResidualThresh(LocalPeakInd(1:nPeaks_Iter));
MaxResidueReduction = MaxResidueReduction(Keep(1:nPeaks_Iter));

LocalPeakInd = LocalPeakInd(Keep);
PeakCodeNo = PeakCodeNo(Keep);
PeakCoef = PeakCoef(Keep,:);
nPeaks = length(PeakCodeNo);
nPeaks_Iter = nPeaks_Iter-sum(~Keep);

%Keep track of which iteration each spot found.
PeakIter = zeros(nPeaks,1);
PeakIter(1:nPeaks_Iter) = iter;
PeakIter(nPeaks_Iter+1:end) = PrevIter;

if o.Graphics==2 && nargout<5
    figure; imagesc(MaxProject);
    colormap(gca,bluewhitered);
    hold on
    scatter(LocalYX(LocalPeakInd(1:nPeaks_Iter),2),LocalYX(LocalPeakInd(1:nPeaks_Iter),1),30,'gx','LineWidth',2);
    hold off
end
end
