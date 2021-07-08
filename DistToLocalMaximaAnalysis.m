ToUse = find(GadGroundTruthLogical(o,'Prob') & o.GadPeakSpots==1 & o.GadPeakColor>250);
nSpots = length(ToUse);
Dist2Maxima = zeros(nSpots,o.nBP,o.nRounds);
for s=1:nSpots
    Dist2Maxima(s,:,:) = iss_view_local_maxima(o,234321,1,'Prob',ToUse(s));
    if mod(s,100)==0
        fprintf([num2str(s),'/',num2str(nSpots),'\n']);
    end
end

nCodes = length(o.GeneNames);
MeanDist = zeros(nSpots,nCodes);
for g=1:nCodes
    GeneChannels = str2double(regexp(cell2mat(o.CharCodes(g)),'\d','match'))+1;    
    CodeIndex = sub2ind([o.nBP,o.nRounds],GeneChannels,1:o.nRounds);
    MeanDist(:,g) = mean(Dist2Maxima(:,CodeIndex),2);
end