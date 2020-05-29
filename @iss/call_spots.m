function o = call_spots(o)
% o = o.call_spots
% calls spots to codes for in-situ sequencing. Run this after find_spots
% 
% produces SpotCodeNo(Spot): gene index for each spot
% SpotScore(Spot): score saying how well the code fits 
% (0...o.nRounds if o.CallSpotsCodeNorm=='Round' else 0...1)
% SpotIntensity(Spot): RMS intensity of the spot
% 
% Using o.UseChannels and o.UseRounds, you can do spot calling
% without using certain rounds and colour channels.
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end
%%

%Only using channels and rounds given by o.UseChannels and o.UseRounds
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end
    
if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);
%o.cSpotColors = o.cSpotColors(:,o.UseChannels,o.UseRounds);

%FILTER OUT REALLY HIGH VALUES
%Good = all(o.cSpotColors(:,:)<10000,2);         
%o.cSpotColors = o.cSpotColors(Good,:,:);
%o.cSpotIsolated = o.cSpotIsolated(Good);
%o.SpotGlobalYX = o.SpotGlobalYX(Good,:);


%Normalise each colour channel by a percentile as to correct for weaker
%colour channels
if strcmpi(o.BleedMatrixType,'Separate')
    p = prctile(o.cSpotColors, o.SpotNormPrctile);
elseif strcmpi(o.BleedMatrixType,'Single')
    p = zeros(1,o.nBP,o.nRounds);
    for b = 1:o.nBP
        bSpotColors = o.cSpotColors(:,b,:);
        p(:,b,:) = prctile(bSpotColors(:), o.SpotNormPrctile);
    end    
else
    warning('Wrong o.BleedMatrixType entry, should be either Separate or Single')
end

pOriginal = p;
HackNo=1;       %1 alters weakest channels, 2 alters most intense
pScale = median(p(:))/10;
DiagMeasure = 0;
nTries = 1;
while DiagMeasure<nChans && nTries<nChans
    SpotColors = bsxfun(@rdivide, o.cSpotColors, p);
    
    % now we cluster the intensity vectors to estimate the Bleed Matrix
    BleedMatrix = zeros(nChans,nChans,nRounds); % (Measured, Real, Round)
    if strcmpi(o.BleedMatrixType,'Separate')
        for r=o.UseRounds
            m = squeeze(SpotColors(o.cSpotIsolated,o.UseChannels,r)); % data: nCodes by nBases
            
            [Cluster, v, s2] = ScaledKMeans(m, eye(nChans));
            for i=1:nChans
                BleedMatrix(:,i,find(o.UseRounds==r)) = v(i,:) * sqrt(s2(i));
            end
        end
        
    elseif strcmpi(o.BleedMatrixType,'Single')
        m = permute(squeeze(squeeze(SpotColors(o.cSpotIsolated,o.UseChannels,o.UseRounds))),[1 3 2]);
        m = squeeze(reshape(m,[],size(m,1)*nRounds,nChans));
        [Cluster, v, s2] = ScaledKMeans(m, eye(nChans));
        for i=1:nChans
            BleedMatrix(:,i,1) = v(i,:) * sqrt(s2(i));
        end
        for r=2:nRounds
            BleedMatrix(:,:,r) = BleedMatrix(:,:,1);
        end
        
    else
        warning('Wrong o.BleedMatrixType entry, should be either Separate or Single')
    end
    
    if o.Graphics
        figure(98043715+nTries); clf
        for i=1:nRounds
            subplot(ceil(nRounds/3),3,i);
            imagesc(BleedMatrix(:,:,i));
            caxis([0 1]);
            title(sprintf('Round %d', o.UseRounds(i)));
            set(gca, 'xtick', 1:nChans);
            set(gca, 'XTickLabel', o.bpLabels(o.UseChannels));
            set(gca, 'ytick', 1:nChans);
            set(gca, 'yTickLabel', o.bpLabels(o.UseChannels));
            if i==4
                xlabel('Actual')
                ylabel('Measured');
            end
        end
        drawnow;
        %     subplot(2,3,6);
        %     caxis([0 1]);
        %     axis off
        %     colormap hot
        % %     colorbar
    end
    
    %If bleed matrix not diagonal, try modifying percentiles of weakest
    %channels
    [~,CurrentBleedMatrixMaxChannel] = max(BleedMatrix(:,:,1));
    DiagMeasure = sum(CurrentBleedMatrixMaxChannel==1:nChans);      %In column i, max square should be in row i if diagonal
    if HackNo==1
        [~,ChangeIntensityChannel] = min(mean(squeeze(p)'));
        p(:,ChangeIntensityChannel,:) = p(:,ChangeIntensityChannel,:)*pScale;
    elseif HackNo==2
        [~,ChangeIntensityChannel] = max(mean(squeeze(p)'));
        p(:,ChangeIntensityChannel,:) = p(:,ChangeIntensityChannel,:)*pScale;
    end
    if DiagMeasure<nChans
        warning('Bleed matrix not diagonal - modifying percentile of channel '+string(ChangeIntensityChannel-1))
    elseif DiagMeasure>=nChans && nTries>1
        fprintf('Bleed matrix now diagonal\n');
    end
    nTries = nTries+1;
    if nTries==4
        p = pOriginal;
        pScale = 0.7;
        HackNo = 2;
        [~,ChangeIntensityChannel] = max(mean(squeeze(p)'));
        p(:,ChangeIntensityChannel,:) = p(:,ChangeIntensityChannel,:)*pScale;
    end
end
if DiagMeasure<nChans
    error('Bleed matrix not diagonal')
end
o.BleedMatrix = BleedMatrix;

% now load in the code book and apply bleeds to it
%codebook_raw = importdata(o.CodeFile);
%CharCode = codebook_raw.textdata(2:end,5);
%GeneName = codebook_raw.textdata(2:end,3);
GeneName = {};
CharCode = {};
fp = fopen(o.CodeFile, 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName=tmp{1};
CharCode=tmp{2};
fclose(fp);

nCodes = size(CharCode,1) - nnz(cellfun(@(v) strcmp(v(1:2),'SW'), CharCode)); % bit of a hack to get rid of Sst and Npy (assume always in the end)

% put them into object o but without the extras
o.CharCodes=CharCode(1:nCodes);
o.GeneNames=GeneName(1:nCodes);

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r=1:o.nRounds
    if r<=o.nRounds-o.nRedundantRounds
        for c=1:nCodes
            [~, NumericalCode(c,r)] = ismember(CharCode{c}(r), o.bpLabels);
        end
    else
        % redundant round - compute codes automatically
        % find pseudobases for this code
        for c=1:nCodes
            PseudoCode = repmat('0',1,o.nRounds-o.nRedundantRounds);
            for p = 1:length(o.RedundantPseudobases)
                PseudoCode(1,ismember(CharCode{c}, o.RedundantPseudobases{p}))=('0'+p);
            end
            % now match them to the redundant codes
            for cc=1:nChans
                rrn = r-o.nRounds+o.nRedundantRounds;
                if ~isempty(regexp(PseudoCode, o.RedundantCodes{rrn,cc}, 'once'))
                    NumericalCode(c,r)=cc;
                end
            end
        end
    end
end

% for r = 1:o.nRounds
% %     if o.AnchorChannel == 2
%         NumericalCode(:,r) = codebook_raw.data(1:nCodes,(r-1)*nChans + (o.AnchorChannel+1:nChans))*(1:o.nBP)';
% %     else
% %         NumericalCode(:,r) = codebook_raw.data(1:nCodes,(r-1)*nChans + (o.DapiChannel+1:nChans-1))*(1:o.nBP)';
% %     end
% end

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:nRounds
        if any(o.UseChannels == NumericalCode(i,o.UseRounds(r))) == 0 continue; end
        BledCodes(i,o.UseChannels+o.nBP*(r-1)) = BleedMatrix(:, find(o.UseChannels == NumericalCode(i,o.UseRounds(r))), r);
        UnbledCodes(i,o.UseChannels(find(o.UseChannels == NumericalCode(i,o.UseRounds(r))))+o.nBP*(r-1)) = 1;
    end
end

%Do code normalisation
FlatSpotColors = SpotColors(:,:);
o.SpotIntensity = sqrt(nansum(FlatSpotColors.^2,2));

if strcmpi(o.CallSpotsCodeNorm,'Round')
    %Normalise so norm in each round is 1
    NormBledCodes = reshape(BledCodes,[nCodes,o.nBP,nRounds]);
    for g=1:nCodes
        for r=1:nRounds
            NormBledCodes(g,:,o.UseRounds(r)) = NormBledCodes(g,:,o.UseRounds(r))/...
                norm(squeeze(NormBledCodes(g,:,o.UseRounds(r))));
        end
    end
    NormBledCodes = reshape(NormBledCodes,[nCodes,o.nBP*nRounds]);
    
    nSpots = size(SpotColors,1); 
    NormSpotColors = SpotColors;
    for s=1:nSpots
        for r=1:nRounds
            NormSpotColors(s,:,o.UseRounds(r)) = NormSpotColors(s,:,o.UseRounds(r))/...
                norm(squeeze(NormSpotColors(s,:,o.UseRounds(r))));
        end
    end
    NormFlatSpotColors = NormSpotColors(:,:);
    
else
    NormBledCodes = bsxfun(@rdivide, BledCodes, sqrt(sum(BledCodes.^2,2)));
    NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, o.SpotIntensity);    
end
%Get rid of NaN values
NormFlatSpotColors(isnan(NormFlatSpotColors)) = 0;
NormBledCodes(isnan(NormBledCodes)) = 0;
SpotScores = NormFlatSpotColors * NormBledCodes';


%Store deviation in spot scores - can rule out matches based on a low
%deviation.
nSpots = size(SpotScores,1);
o.SpotScoreDev = zeros(nSpots,1);
for s=1:nSpots
    o.SpotScoreDev(s) = std(SpotScores(s,:));
end

[o.SpotScore, BestCode] = max(SpotScores,[],2);
o.SpotCodeNo = uint16(BestCode);
o.SpotCombi = true(size(o.SpotCodeNo,1),1);

o.BledCodes = BledCodes;
o.UnbledCodes = UnbledCodes;
o.NormBledCodes = NormBledCodes;
o.cNormSpotColors = reshape(NormFlatSpotColors,size(o.cSpotColors));
