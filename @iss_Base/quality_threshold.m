function QualOK = quality_threshold(o,Method)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold
%Method = 'DotProduct','Prob' or 'Pixel' to consider gene assignments given
%by o.SpotCodeNo, o.pSpotCodeNo and o.pxSpotCodeNo respectively.

if nargin<2 || isempty(Method)
    Method = 'DotProduct';
end

if ~ismember({Method},o.CallMethods)
    error('Method invalid, must be member of o.CallMethods.');
end
pf = o.CallMethodPrefix(Method);

if strcmpi('Prob',Method) || strcmpi('Pixel',Method) || strcmpi('GroundTruth',Method)
%     QualOK = (o.([pf,'SpotScore'])>o.pScoreThresh & o.([pf,'SpotIntensity'])>o.pIntensityThresh2 | ...
%     o.([pf,'SpotScore'])>o.pScoreThresh2 & o.([pf,'SpotScore'])+o.([pf,'LogProbOverBackground'])>o.pLogProbThresh2 &...
%     o.([pf,'SpotIntensity'])>o.pIntensityThresh); 
%     QualOK = QualOK & o.([pf,'SpotScore'])>=0;       %3rd best spots do more harm than good.
    %QualOK = QualOK & o.pSpotIntensity2 > o.pIntensity2Thresh;
    QualOK = o.([pf,'SpotScore'])>0 & o.([pf,'LogProbOverBackground'])+o.pQualParam1*o.([pf,'SpotScore'])>o.pQualThresh1 | ...
       o.([pf,'SpotScore'])==0 & o.([pf,'LogProbOverBackground'])+o.pQualParam2*o.([pf,'SpotScore'])>o.pQualThresh2 | ...
       o.([pf,'SpotScore'])<0 & o.([pf,'LogProbOverBackground'])>o.pQualThresh3 & o.([pf,'SpotScore'])>o.pQualThresh4;
    QualOK = QualOK & o.([pf,'SpotIntensity'])>o.pIntensityThresh;
elseif strcmpi('DotProduct',Method)
    QualOK = o.([pf,'SpotCombi']) & o.([pf,'SpotScore'])>o.CombiQualThresh &...
        o.([pf,'SpotIntensity'])>o.CombiIntensityThresh & o.([pf,'SpotScoreDev'])>o.CombiDevThresh;
elseif strcmpi('OMP',Method)
     %Old method below
     %QualOK = o.ompNeighbNonZeros>o.ompNeighbThresh | (o.ompSpotIntensity>o.ompIntensityThresh & o.ompNeighbNonZeros>o.ompNeighbThresh2);
     %QualOK = QualOK & o.ompSpotIntensity2 > o.ompIntensity2Thresh;
     %New method, found using PyTorch
     QualOK = o.([pf,'NeighbNonZeros'])>o.ompNeighbThresh | o.([pf,'SpotIntensity'])>o.ompIntensityThresh |...
         o.([pf,'SpotScore'])>o.ompScoreThresh;
elseif strcmpi('Spatial',Method)
    QualOK = o.([pf,'SpotScore'])>0;        %All spots at the moment.
else
    %If new method, just accept everything
    QualOK = o.([pf,'SpotCodeNo'])>0;   
end


% % HACK ALERT
% QualOK = QualOK & o.cSpotIsolated;

% nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));
% 
% % now extras - they have their own thresholds, set manually for each type
% for i=1:size(o.ExtraCodes,1)
%     MySpots = (o.SpotCodeNo == nCombiCodes+i);
%     QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
% end