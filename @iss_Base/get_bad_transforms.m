function [A, D, D_average, PcFailed] = get_bad_transforms(o,NonemptyTiles,nMatches,D)
%% [A, D, D_average, PcFailed] = get_bad_transforms(o,NonemptyTiles,nMatches,D)
%
% This finds the transforms that need to be re-evaluated. It also updates
% the starting transforms for this re-evaluation to be the average
% transform. The regularization will also ensure the new transform found is
% close to this average transform.
%
% o: iss object
% NonemptyTiles: tile numbers of tiles to evaluate
% nMatches(t,b,r): the number of matches found with the current transform
% for tile t, channel b, round r.
% D(:,:,t,r,b): the current transform found for tile t, channel b, round r.
% A(b,:): the median chromatic aberration scaling for channel b.
% D: the current transforms with the ones for which PcFailed set to D_average.
% D_average: the average transform - scaling is A(b,:). Shift is average
% over all colour channels for that tile and round. Rotation=0.
% PcFailed(t,b,r): This is true if tile t, channel b, round r is to be
% re-evaluated.

%% Find bad transforms that need to be re-evaluated using regularization
%These either have few matches, bad CA scalings or bad shifts
%Require or condition below as sometimes one spot can match to 1000s if
%scaling is set to 0.
BadMatches = nMatches<o.PcMinSpots | o.AllBaseSpotNo<o.PcMinSpots;
PcFailed = BadMatches;

CA_Scale1 = permute(squeeze(D(1,1,:,:,:)),[1,3,2]);
CA_Scale2 = permute(squeeze(D(2,2,:,:,:)),[1,3,2]);
Scale_median1 = zeros(size(BadMatches));
Scale_median2 = zeros(size(BadMatches));
for b=o.UseChannels
    CA_ToUse = false(size(BadMatches));
    CA_ToUse(:,b,:)=true;
    Scale_median1(:,b,:) = median(squeeze(CA_Scale1(~PcFailed&CA_ToUse)));
    Scale_median2(:,b,:) = median(squeeze(CA_Scale2(~PcFailed&CA_ToUse)));
end
BadCA = abs(CA_Scale1-Scale_median1)>o.PcMaxScaleDev | ...
    abs(CA_Scale2-Scale_median2)>o.PcMaxScaleDev;
BadCA = BadCA & nMatches<o.PcMinSpotsScaling;

    
PcFailed = BadMatches | BadCA;
Shift_median1 = zeros(size(BadMatches));
Shift_median2 = zeros(size(BadMatches));
for t=NonemptyTiles
    for r = o.UseRounds
        Shift_t = median(squeeze(D(3,:,t,r,~PcFailed(t,:,r)))');
        Shift_median1(t,:,r) = Shift_t(1);
        Shift_median2(t,:,r) = Shift_t(2);
    end
end
BadShift = abs(permute(squeeze(D(3,1,:,:,:)),[1,3,2])-Shift_median1)>o.PcMaxShiftDev |...
           abs(permute(squeeze(D(3,2,:,:,:)),[1,3,2])-Shift_median2)>o.PcMaxShiftDev;
       
PcFailed = BadMatches | BadCA | BadShift;

%% Update starting transforms of these transforms
%Record chromatic aberration by taking median of all transforms that worked
A = zeros(o.nBP,2);
for b=o.UseChannels
    CA_ToUse = false(size(BadMatches));
    CA_ToUse(:,b,:)=true;
    A(b,1) = median(squeeze(CA_Scale1(~PcFailed&CA_ToUse)));
    A(b,2) = median(squeeze(CA_Scale2(~PcFailed&CA_ToUse)));
end

D_average = zeros(size(D));
for t=NonemptyTiles
    for r = o.UseRounds
        Shift_t = median(squeeze(D(3,:,t,r,~PcFailed(t,:,r)))');
        for b=o.UseChannels
            D_average(1,1,t,r,b) = A(b,1);
            D_average(2,2,t,r,b) = A(b,2);
            D_average(3,:,t,r,b) = Shift_t;
            if PcFailed(t,b,r)
                %Set starting transform to average transforms for those
                %that are to be re-evaluated.
                D(:,:,t,r,b) = D_average(:,:,t,r,b);
            end
        end
    end
end
