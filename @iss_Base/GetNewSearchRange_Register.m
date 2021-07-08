function o = GetNewSearchRange_Register(o,shifts,direction)

%For a particular direction, this considers all shifts found so far
%and adjusts the search range based on them.

%Find extremes of shifts found so far
MinY = min(shifts(:,1));
MaxY = max(shifts(:,1));
MinX = min(shifts(:,2));
MaxX = max(shifts(:,2));

%Adjust search ranges
if strcmpi(direction,'South')
    o.RegSearch.South.Y = MinY-2*o.RegStep(1):o.RegStep(1):MaxY+2*o.RegStep(1);
    o.RegSearch.South.X = MinX-2*o.RegStep(2):o.RegStep(2):MaxX+2*o.RegStep(2);
elseif strcmpi(direction,'East')
    o.RegSearch.East.Y = MinY-2*o.RegStep(1):o.RegStep(1):MaxY+2*o.RegStep(1);
    o.RegSearch.East.X = MinX-2*o.RegStep(2):o.RegStep(2):MaxX+2*o.RegStep(2);
end
end
