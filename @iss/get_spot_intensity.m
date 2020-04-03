function SpotIntensity = get_spot_intensity(o,SpotCodeNo)
%This gives a modified spot intensity taking account of the gene it is
%assigned to.
%For spot s, it takes all colour channels that appear in the gene
%SpotCodeNo(s). Call this SpotCode2, then SpotIntensity(s) =
%mean(SpotCode2(InCharCode))-mean(SpotCode2(NotInCharCode))
%Hence a high SpotIntensity indicates a high intensity and a good match.

nCodes = length(o.GeneNames);
CodeIndex = zeros(nCodes,o.nRounds);
NonCodeIndex = cell(nCodes,1);
for g=1:nCodes
    GeneChannels = str2double(regexp(cell2mat(o.CharCodes(g)),'\d','match'))+1;    
    CodeIndex(g,:) = sub2ind([o.nBP,o.nRounds],GeneChannels,1:o.nRounds);
    UnusedChannels = setdiff(1:o.nBP,GeneChannels);
    UnusedChannelIndex = sub2ind([o.nBP,o.nRounds],repelem(UnusedChannels,1,o.nRounds),repmat(1:o.nRounds,1,length(UnusedChannels)));
    NonCodeIndex{g} = setdiff(1:o.nRounds*o.nBP,[CodeIndex(g,:),UnusedChannelIndex]);
end

nSpots = length(SpotCodeNo);
SpotIntensity = zeros(nSpots,1);

fprintf('Percentage of spot intensities found:       ');
for s=1:nSpots
    SpotCode = o.cSpotColors(s,:);
    SpotIntensity(s) = mean(SpotCode(CodeIndex(SpotCodeNo(s),:)))-mean(SpotCode(NonCodeIndex{SpotCodeNo(s)}));
    if mod(s,round(nSpots/100))==0
        Percent = sprintf('%.6f', round(s*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end
fprintf('\n');
end


% %This was original method used, differs slightly as when colour channel
% %appeared in more than one round, took mean of all rounds it appeared first. 
% %This method was much slower though and pretty much same results
% nSpots = length(SpotCodeNo);
% SpotIntensity = zeros(nSpots,1);
% RoundCode = 1:o.nRounds;
% for s=1:nSpots
%     SpotCode = o.cSpotColors(s,:,:);
%     numCharCode = str2double(regexp(cell2mat(o.CharCodes(SpotCodeNo(s))),'\d','match'))+1;    
%     sIntensity = zeros(7,1);
%     for b=1:o.nBP
%         UseRounds = find(numCharCode==b);
%         if ~isempty(UseRounds)
%             sIntensity(b)= mean(SpotCode(:,b,UseRounds))-mean(SpotCode(:,b,setdiff(RoundCode,UseRounds)));
%         end
%     end
%     SpotIntensity(s) = mean(sIntensity(sIntensity~=0));
% end
% 
