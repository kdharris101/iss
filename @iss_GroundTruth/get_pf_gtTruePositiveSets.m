function o = get_pf_gtTruePositiveSets(o,Method)
%% o = o.get_pf_gtTruePositiveSets(Method)
% True positive if spot found by Method assigned as ground truth gene and 
% within o.gtTruePosMaxSep of ground truth round/channel peak.
% False positive if spot found by Method as ground truth gene but further 
% from ground truth round/channel peak than o.gtFalsePosMinSep. 
% and has ground truth round/channel intensity below o.gtColorFalsePositiveThresh

if ~ismember({Method},setdiff(o.CallMethods,{'DotProduct'}))
    error('Method invalid, must be member of o.CallMethods (Not DotProduct).');
end
pf = o.CallMethodPrefix(Method);
SetSize = size(o.([pf,'SpotCodeNo']));
o.([pf,'_gtIdentity']) = cell(o.nRounds+o.nExtraRounds,o.nBP);
o.([pf,'_gtFound']) = cell(o.nRounds+o.nExtraRounds,o.nBP);

%get peaks in the ground truth rounds/channels
o.gtTruePositiveSet = get_gtTruePositiveSet(o,o.gtColorTruePositiveThresh,true);
%For FalsePositiveSet, relax the constraint that peaks are only high in
%one channel and also use lower threshold.
gtFalsePositiveSet = o.get_gtTruePositiveSet(o.gtColorFalsePositiveThresh,false);

for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        gtSpotGlobalYX_TP = o.gtSpotGlobalYX{r,b}(o.gtTruePositiveSet{r,b},:);
        fprintf('There are %d %s peak spots\n', length(gtSpotGlobalYX_TP),...
            o.GeneNames{o.gtGeneNo(r,b)});
        
        gtSpotGlobalYX_FP = o.gtSpotGlobalYX{r,b}(gtFalsePositiveSet{r,b},:);
        
        %For each point GadPeakTruePositive, find nearest Gad spot found by pixel
        %based method
        pfGeneSpotsIndex = find(o.([pf,'SpotCodeNo'])==o.gtGeneNo(r,b));
        pfGeneSpotsYX = o.([pf,'SpotGlobalYX'])(pfGeneSpotsIndex,:);
        tree_TP = KDTreeSearcher(pfGeneSpotsYX);
        [Index2_TP,Dist_TP] = tree_TP.knnsearch(gtSpotGlobalYX_TP);
        
        %Spot is true positive if within thresh disance of a peak in true
        %postive round/channel.
        %I.e. of the 2980 true positives peaks in ground truth channel/round, 
        %how many spots identified as ground truth gene using pf gene calling method
        %are close to these.
        TruePositiveIndex = unique(pfGeneSpotsIndex(Index2_TP(Dist_TP<o.gtTruePosMaxSep)));
        
        %False positive if spot is assigned as ground truth gene using pf gene calling 
        %method but is further from ground truth peak than thresh distance and has 
        %ground truth channel/round intensity below threshold.
        treeFalsePositive = KDTreeSearcher(gtSpotGlobalYX_FP);
        [Index2_FP,Dist_FP] = treeFalsePositive.knnsearch(pfGeneSpotsYX);
        FalsePositiveIndex = unique(pfGeneSpotsIndex(Dist_FP>o.gtFalsePosMinSep));
        
        o.([pf,'_gtIdentity']){r,b} = zeros(SetSize);
        o.([pf,'_gtIdentity']){r,b}(TruePositiveIndex) = 1;
        FalsePosSet = false(SetSize);
        FalsePosSet(FalsePositiveIndex) = true;
        FalsePosSet = FalsePosSet &...
            o.([pf,'_gtColor'])(:,b,r)<o.gtColorFalsePositiveThresh;
        o.([pf,'_gtIdentity']){r,b}(FalsePosSet) = 2;
        
        %Save ground truth that were found (1) and those missed (2). 
        FOUNDTruePositive = Dist_TP<o.gtTruePosMaxSep;
        MISSEDTruePositive = Dist_TP>=o.gtTruePosMaxSep;
        TruePositiveSetIndex = find(o.gtTruePositiveSet{r,b});
        o.([pf,'_gtFound']){r,b} = zeros(size(o.gtTruePositiveSet{b}));
        o.([pf,'_gtFound']){r,b}(TruePositiveSetIndex(FOUNDTruePositive))=1;
        o.([pf,'_gtFound']){r,b}(TruePositiveSetIndex(MISSEDTruePositive))=2;
    end
end

end