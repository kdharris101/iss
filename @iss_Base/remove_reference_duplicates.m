function o = remove_reference_duplicates(o,NonemptyTiles)
%Reference spots are made up of spots on multiple colour channels so you
%can get duplicates. This removes these by only keeping spots in a colour
%channel if closest spot in all other colour channels is further away than
%o.DetectionRadius
%Note, this preferentially keeps spots on latter colour channels. Shouldn't
%matter too much if chromatic aberration is small though. In future, maybe
%keep spots on colour channel in which intensity is greatest.

for t=NonemptyTiles
    for b = setdiff(o.ReferenceSpotChannels,o.ReferenceChannel)
        AnchorSpots = vertcat(o.RawLocalYX{t,setdiff(o.ReferenceSpotChannels,b)});
        k0 = KDTreeSearcher(AnchorSpots);
        [~, dist] = k0.knnsearch(o.RawLocalYX{t,b}, 'k', 1);
        o.RawLocalYX(t,b) = {o.RawLocalYX{t,b}(dist>o.DetectionRadius,:)};
        o.RawIsolated(t,b) = {o.RawIsolated{t,b}(dist>o.DetectionRadius,:)};
    end
end

end