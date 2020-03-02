function BledCodes = change_bled_codes(o,NewBleedMatrix)
%Given a new bleed matrix, this returns the corresponding bled codes

nCodes = length(o.GeneNames);
% create numerical code (e.g. 33244 for CCGAA)
FullCharCodes = cell2mat(o.CharCodes)';
numCharCode = str2double(regexp(FullCharCodes(:)','\d','match'))+1;
NumericalCode = reshape(numCharCode,[o.nRounds,nCodes])';

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:o.nRounds
        if any(o.UseChannels == NumericalCode(i,o.UseRounds(r))) == 0 continue; end
        BledCodes(i,o.UseChannels+o.nBP*(r-1)) = NewBleedMatrix(:, find(o.UseChannels == NumericalCode(i,o.UseRounds(r))), r);
    end
end

% if o.ArtificialGenes
%     BledCodes(74:end,:) = o.pBledCodes(74:end,:); %Artificial genes remain same
% end