function o = AmendPCR(o)
%% o = AmendPCR(o)
%
% This function takes all tiles/rounds/channels for which
% o.nMatches(t,b,r)<o.MinSpots as indicated by o.PcFailed. It then assumes
% that the transform found for these are erronous because it will have
% overfitted to the few matches. Instead, it finds a new transform by
% assuming the mean chromatic aberration o.A(b) and a shift found using a
% Fft method on the images. The rotation is set to 0. 
%
% Alters o.D. Previous transforms are now saved as o.FindSpotsInfo.DFailed

%%
NonemptyTiles = find(~o.EmptyTiles)';
if ~isfield(o.FindSpotsInfo,'DFailed')
    o.FindSpotsInfo.DFailed = o.D;
end

fprintf('PCR - Correcting transforms for tile   ');
for t=NonemptyTiles
    if t<10
        fprintf('\b%d',t);
    else
        fprintf('\b\b%d',t);
    end
    for b=o.UseChannels
        for r=o.UseRounds
            if ~o.PcFailed(t,b,r); continue; end
            %Use mean chromatic aberration for colour channel and forget
            %rotation
            o.D(:,:,t,r,b) = o.D(:,:,t,r,b)*0;
            o.D(1,1,t,r,b) = o.A(b,1);
            o.D(2,2,t,r,b) = o.A(b,2);
            
            %Obtain shift through Fft method. Only use this if it is
            %comparable to shifts found for other channels in round, for
            %which PCR worked. Otherwise set to shift to colour channel
            %that has chromatic aberration that is closest.
            RoundShift = squeeze(o.D(3,:,t,r,~o.PcFailed(t,:,r)))';   
            MinAllowedShift = min(RoundShift)-std(RoundShift);
            MaxAllowedShift = max(RoundShift)+std(RoundShift);
            MinAllowedShift = -Inf;
            MaxAllowedShift = Inf;
            [shift, ~] = o.get_Fft_shift_single(t,r,b,...
                t,o.ReferenceRound,o.ReferenceChannel,'SubPixel');
            if all(shift<MaxAllowedShift) && all(shift>MinAllowedShift)
                o.D(3,:,t,r,b) = shift;
            else
                [~,BestChannel] = min(sum(abs(o.A./~o.PcFailed(t,:,r)'-o.A(b,:)),2));
                o.D(3,:,t,r,b) = o.D(3,:,t,r,BestChannel);
                warning(['Tile %d, channel %d, round %d\nShift found: [%.2f,%.2f]\n'...
                    'This is outside the range between [%.2f, %.2f] and [%.2f, %.2f]\n'...
                    'Instead using shift of colour channel %d: [%.2f,%.2f]'],...
                    t,b,r,shift(1),shift(2),MinAllowedShift(1),MinAllowedShift(2),...
                    MaxAllowedShift(1), MaxAllowedShift(2),BestChannel,o.D(3,1,t,r,b),o.D(3,2,t,r,b));
            end            
            
        end
    end
end
fprintf('\n');
