function o = call_spots(o)
% o = o.call_spots
% calls spots to codes for in-situ sequencing. Run this after find_spots
% 
% produces SpotGene{Spot}: name of gene for each spot
% SpotCode{Spot}: text representation of code for each spot 
% SpotScore(Spot): score saying how well the code fits (0...1)
% SpotIntensity(Spot): RMS intensity of the spot
% 
% Using o.UseChannels and o.UseRounds, you can do spot calling
% without using certain rounds and colour channels.
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

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

%Filter out high values in problematic round 3 colour 1 in bottom right
%Bad = o.SpotGlobalYX(:,1) > 479 & o.SpotGlobalYX(:,1) < 1878 & o.SpotGlobalYX(:,2)>6695 &...
%    o.cSpotColors(:,1,3)>1000;      
%o.cSpotColors = o.cSpotColors(Bad==0,:,:);
%o.cSpotIsolated = o.cSpotIsolated(Bad==0);
%o.SpotGlobalYX = o.SpotGlobalYX(Bad==0,:);

SpotColors = bsxfun(@rdivide, o.cSpotColors, prctile(o.cSpotColors, o.SpotNormPrctile));

% now we cluster the intensity vectors to estimate the Bleed Matrix
BleedMatrix = zeros(nChans,nChans,nRounds); % (Measured, Real, Round)
for r =o.UseRounds
    m = squeeze(SpotColors(o.cSpotIsolated,o.UseChannels,r)); % data: nCodes by nBases
    
    [Cluster, v, s2] = ScaledKMeans(m, eye(nChans));
    for i=1:nChans
        BleedMatrix(:,i,r) = v(i,:) * sqrt(s2(i));
    end
end

if o.Graphics
    figure(98043765); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i); 
        imagesc(BleedMatrix(:,:,i)); 
        caxis([0 1]); 
        title(sprintf('Cycle %d', o.UseRounds(i))); 
        set(gca, 'xtick', 1:nChans);
        set(gca, 'XTickLabel', o.bpLabels(o.UseChannels));
        set(gca, 'ytick', 1:nChans);
        set(gca, 'yTickLabel', o.bpLabels(o.UseChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
%     subplot(2,3,6);
%     caxis([0 1]); 
%     axis off
%     colormap hot
% %     colorbar
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

if 1 % 0 to just use original codes
    NormBledCodes = bsxfun(@rdivide, BledCodes, sqrt(sum(BledCodes.^2,2)));
    FlatSpotColors = SpotColors(:,:);
    o.SpotIntensity = sqrt(nansum(FlatSpotColors.^2,2));
    NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, o.SpotIntensity);
    
    %Get rid of NaN values
    NormFlatSpotColors(isnan(NormFlatSpotColors)) = 0;
    NormBledCodes(isnan(NormBledCodes)) = 0;    
    SpotScores = NormFlatSpotColors * NormBledCodes';
    
else
    % HACK ALERT
    NormBledCodes = bsxfun(@rdivide, BledCodes(:,1:20), sqrt(sum(BledCodes(:,1:20).^2,2)));
    FlatSpotColors = SpotColors(:,1:20);
    o.SpotIntensity = sqrt(sum(FlatSpotColors.^2,2));
    NormFlatSpotColors = bsxfun(@rdivide, FlatSpotColors, o.SpotIntensity);

    SpotScores = NormFlatSpotColors * NormBledCodes';
end

[o.SpotScore, BestCode] = max(SpotScores,[],2);
o.SpotCodeNo = uint16(BestCode);
o.SpotCombi = true(size(o.SpotCodeNo,1),1);

o.BledCodes = BledCodes;
o.UnbledCodes = UnbledCodes;
o.NormBledCodes = NormBledCodes;
o.cNormSpotColors = reshape(NormFlatSpotColors,size(o.cSpotColors));
