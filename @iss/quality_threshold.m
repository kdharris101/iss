function QualOK = quality_threshold(o)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold

QualOK = (o.SpotCombi & o.SpotScore>o.CombiQualThresh & o.SpotIntensity>o.CombiIntensityThresh);

AnchorsOK = ones(size(QualOK));
AnchorsOK(o.SpotCombi) = (sum(o.cAnchorIntensities>o.DetectionThresh,2)>o.CombiAnchorsReq);

QualOK = QualOK & AnchorsOK;

nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));

% now extras - they have their own thresholds, set manually for each type
for i=1:size(o.ExtraCodes,1)
    MySpots = (o.SpotCodeNo == nCombiCodes+i);
    QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
end