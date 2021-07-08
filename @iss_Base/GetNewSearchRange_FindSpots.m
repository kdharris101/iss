function o = GetNewSearchRange_FindSpots(o,t,r)

%For a particular round r, this considers all shifts found so far up to
%tile t and adjusts the search range based on them.

%Find extremes of shifts found so far
MinY = min(o.D0(1:t,1,r));
MaxY = max(o.D0(1:t,1,r));
MinX = min(o.D0(1:t,2,r));
MaxX = max(o.D0(1:t,2,r));

%Adjust search ranges
o.FindSpotsSearch{r}.Y = MinY-2*o.FindSpotsStep(1):o.FindSpotsStep(1):MaxY+2*o.FindSpotsStep(1);
o.FindSpotsSearch{r}.X = MinX-2*o.FindSpotsStep(2):o.FindSpotsStep(2):MaxX+2*o.FindSpotsStep(2);

end