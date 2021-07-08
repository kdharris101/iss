function [AllSpotCodeNo,AllSpotCoef,AllSpotLocalInd,AllSpotIter,AllResidueReduction,...
    IterSpotColors,IterConvImages] = get_spatial_spots(o,SpotColors,LocalYX)
%% [AllSpotCodeNo,AllSpotCoef,AllSpotLocalInd,AllSpotIter,AllResidueReduction,...
%  IterSpotColors,IterConvImages] = o.get_spatial_spots(SpotColors,LocalYX)
%  Finds peaks in convolution images and returns info about peaks.
%  SpotColors: z-scored spot colors.
%  LocalYX: yx coordinate of these spots on tile.
%  AllSpotCodeNo: gene found at LocalYX(AllSpotLocalInd,:) on tile.
%  AllSpotCoef(s,a) is the coefficient for spatial eigenvector a for spot s.
%  AllSpotIter(s): iteration that spot s was found on.
%  AllResidueReduction(s): reduction in residue caused by the inclusion of
%  gene AllSpotCodeNo(s) at LocalYX(AllSpotLocalInd(s),:).
%  IterSpotColors(it+1) is SpotColors after it iterations.
%  IterSpotColors(1)=SpotColors.
%  IterConvImages{it}(:,:,g) is the convolution image for gene g with the
%  first spatial eigenvector at iteration it. 
%% 
if nargout>5
    IterSpotColors = cell(o.sp_nIter+1,1);
    IterSpotColors{1} = SpotColors;
end
if nargout>6
    IterConvImages = cell(o.sp_nIter,1);
end
nCodes = length(o.CharCodes);
BackgroundCodeNo = nCodes+1:size(o.spBledCodes,1);
AllSpotCodeNo = cell(o.sp_nIter,1);
AllSpotCoef = cell(o.sp_nIter,1);
AllSpotLocalInd = cell(o.sp_nIter,1);
AllSpotIter = cell(o.sp_nIter,1);
AllResidueReduction = cell(o.sp_nIter,1);

%OMP stops when reduction in residual drops below ResidualThresh.
%Prctile bit gets 2nd largest intensity for each spot.
ResidualThresh = prctile(abs(SpotColors(:,:))',47.5*100/49.0)';
ResidualThresh = o.ResidualThreshParam*ResidualThresh;
ResidualThresh(ResidualThresh<o.ResidualThreshMin) = o.ResidualThreshMin;
ResidualThresh(ResidualThresh>o.ResidualThreshMax) = o.ResidualThreshMax;
%ResidualThresh = ResidualThresh+o.ResidualReductionThresh;

fprintf('Iteration    ');
for it=1:o.sp_nIter
    fprintf('\b\b\b%02d ',it);
    %Fit background at each pixel.
    SpotColors = o.update_spot_background(SpotColors,BackgroundCodeNo);
    %Find location of spots on this iteration and update coefficients of
    %spots found on previous iteration. 
    if nargout>6
        [itLocalInd,itSpotCodeNo,itSpotCoef,itSpotIter,AllResidueReduction{it},IterConvImages{it}] =...
            o.get_conv_peaks(LocalYX,SpotColors,ResidualThresh,...
            cell2mat(AllSpotLocalInd(1:it)),cell2mat(AllSpotCodeNo(1:it)),...
            cell2mat(AllSpotIter),it,BackgroundCodeNo);
    else
        [itLocalInd,itSpotCodeNo,itSpotCoef,itSpotIter,AllResidueReduction{it}] =...
            o.get_conv_peaks(LocalYX,SpotColors,ResidualThresh,...
            cell2mat(AllSpotLocalInd(1:it)),cell2mat(AllSpotCodeNo(1:it)),...
            cell2mat(AllSpotIter),it,BackgroundCodeNo);
    end
    %Based on these spots, update spot color at all pixels. 
    SpotColors = o.update_spot_colors(...
        LocalYX,SpotColors,itLocalInd,itSpotCodeNo,itSpotCoef);
    if nargout>5
        IterSpotColors{it+1} = SpotColors;
    end
    AllSpotCodeNo{it} = itSpotCodeNo(itSpotIter==it);
    AllSpotLocalInd{it} = itLocalInd(itSpotIter==it);
    AllSpotCoef{it} = itSpotCoef(itSpotIter==it,:);
    if it>1
        for it2=1:it-1
            AllSpotCoef{it2} = AllSpotCoef{it2}+itSpotCoef(itSpotIter==it2,:);
        end
    end
    AllSpotIter{it} = ones(size(AllSpotCodeNo{it}))*it;
end
fprintf('\n');

AllSpotCodeNo = cell2mat(AllSpotCodeNo);
AllSpotCoef = cell2mat(AllSpotCoef);
AllSpotLocalInd =  cell2mat(AllSpotLocalInd);
AllSpotIter = cell2mat(AllSpotIter);
AllResidueReduction = cell2mat(AllResidueReduction);

end

