function o = convert_iss_data(oOld)
%% objNew = convert_iss_data(objOld);
% Converts data from iss_PixelBased branch on GitHub.
% oOld: iss object from issPixelBased branch on GitHub.
% o: iss object from iss_PixelBased class of current iss software.
% Note that you can still use the same LookupTable as made with old class.

%% Update all properties that exist for both classes.
DesiredClass = iss_PixelBased;
o = DesiredClass;
NewClassProp = metaclass(DesiredClass).PropertyList;
for k = 1:length(NewClassProp)
    try
        o.(NewClassProp(k).Name) = oOld.(NewClassProp(k).Name);
    catch
        ...
    end
end

%% Deal with differences in properties/property names
% Extract variables
o.StripHack = false;           %Did not account for strips of zeros in raw data 
                               %with old method.
%% FindSpots output
%Change form of transform
if ~isempty(o.D) && ~isempty(o.A)
   o.PointCloudMethod = 2;        %Used PCR2 in this branch.
   o.FindSpotsInfo.D_fromPCR2 = o.D;
   D = zeros(3,2,size(o.TileOrigin,1),o.nRounds,o.nBP);
   for t=1:size(o.TileOrigin,1)
       for r=1:o.nRounds
           for b=1:o.nBP
               D(:,:,t,r,b) = o.A(b)*o.D(:,:,t,r);
           end
       end
   end
   o.D = D;
end

%% DotProduct outputs
o.dpNormSpotColors = oOld.cNormSpotColors;
o.dpSpotCodeNo = oOld.SpotCodeNo;
o.dpSpotColors = oOld.cSpotColors;
o.dpSpotCombi = oOld.SpotCombi;
o.dpSpotGlobalYX = oOld.SpotGlobalYX;
o.dpSpotIntensity = oOld.SpotIntensity;
o.dpSpotIsolated = oOld.cSpotIsolated;
o.dpSpotScore = oOld.SpotScore;
o.dpSpotScoreDev = oOld.SpotScoreDev;

if ~isempty(o.dpSpotColors)       
    %Normalisation to make Bleed Matrix
    if strcmpi(o.BleedMatrixType,'Separate')
        o.BledCodesPercentile = prctile(o.dpSpotColors, o.SpotNormPrctile);
    elseif strcmpi(o.BleedMatrixType,'Single')
        o.BledCodesPercentile = zeros(1,o.nBP,o.nRounds);
        for b = 1:o.nBP
            bSpotColors = o.dpSpotColors(:,b,:);
            o.BledCodesPercentile(:,b,:) = prctile(bSpotColors(:), o.SpotNormPrctile);
        end
    end   
    o.dpLocalTile = o.get_local_tile(o.dpSpotGlobalYX);
end

%% Prob method outputs
o.pSpotColors = oOld.cSpotColors;
o.pSpotGlobalYX = oOld.SpotGlobalYX;

if ~isempty(o.pSpotCodeNo)
    o.ProbMethod = 2;              %Used old prob method in old branch.
    o.ScoreScale = 1;              %ProbMethod=2 requires ScoreScale=1.
    
    x = min(o.pSpotColors(:))-1:max(o.pSpotColors(:))-1;
    %Get background distribution by shifting Histogram so aligns with 0 of spot colors.
    BackgroundGamma = dirac(x);
    BackgroundGamma(BackgroundGamma==inf)=1;
    o.BackgroundProb = zeros(length(x),o.nBP,o.nRounds);
    for b=1:o.nBP
        for r=1:o.nRounds
            o.BackgroundProb(:,b,r) = conv(BackgroundGamma,o.HistProbs(:,b,r),'same');
        end
    end
    o.pLocalTile = o.get_local_tile(o.pSpotGlobalYX);
end

%% Pixel method ouputs
if ~isempty(o.pxSpotGlobalYX)
    o.pxLocalTile = o.get_local_tile(o.pxSpotGlobalYX);
end

end

