function QualOK = quality_threshold(o,Method)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold
%Method = 'DotProduct','Prob' or 'Pixel' to consider gene assignments given
%by o.SpotCodeNo, o.pSpotCodeNo and o.pxSpotCodeNo respectively.


if strcmpi('Prob',Method)
    QualOK = (o.pSpotScore>o.pScoreThresh & o.pSpotIntensity>0 | ...
    o.pSpotIntensity>o.pIntensityThresh & o.pLogProbOverBackground>o.pLogProbThresh & o.pSpotScore+o.pSpotScoreDev>o.pDevThresh...
    & o.pSpotScore>0); 
%| o.pSpotIntensity>1000);
elseif strcmpi('Pixel',Method)
    QualOK = (o.pxSpotScore>o.pScoreThresh & o.pxSpotIntensity>0 | ...
    o.pxSpotIntensity>o.pIntensityThresh & o.pxLogProbOverBackground>o.pLogProbThresh & o.pxSpotScore+o.pxSpotScoreDev>o.pDevThresh...
    & o.pxSpotScore>0);
%| o.pSpotIntensity>1000);
elseif strcmpi('DotProduct',Method)
    QualOK = (o.SpotCombi & o.SpotScore>o.CombiQualThresh & o.SpotIntensity>o.CombiIntensityThresh & o.SpotScoreDev>o.CombiDevThresh);
else
    error('Method not valid, must be DotProduct, Prob or Pixel');
end


% % HACK ALERT
% QualOK = QualOK & o.cSpotIsolated;

nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));

% now extras - they have their own thresholds, set manually for each type
for i=1:size(o.ExtraCodes,1)
    MySpots = (o.SpotCodeNo == nCombiCodes+i);
    QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
end