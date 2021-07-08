function NewHistCounts = SmoothHistCounts(o)
%Fit where counts fall below PosStartGuess with an exponential decay
%Idea is to avoid large variation in LogProbOverBackground at large
%intensities
StartGuess = [o.SmoothHistCountLimit, o.SmoothHistCountParamGuess];
nBins = size(o.HistCounts,1);
NewHistCounts = o.HistCounts;
for r=1:o.nRounds
    for b=1:o.nBP
        y = o.HistCounts(:,b,r);
        %Find indices closest to centre where count falls below 250.
        Idx = find(y<o.SmoothHistCountLimit);
        ValPos = min(Idx(Idx>o.TilePixelValueShift));
        ValNeg = max(Idx(Idx<o.TilePixelValueShift));
        xPos = ((ValPos:nBins)-ValPos)';        %Start at 0 so initial Pred easier
        fPos = fit(xPos,y(xPos+ValPos),'exp1','StartPoint',StartGuess);
        xNeg = ((1:ValNeg)-1)';                 %Start at 0 so initial Pred easier
        fNeg = fit(xNeg,flip(y(xNeg+1)),'exp1','StartPoint',StartGuess);
        NewHistCounts(xPos+ValPos,b,r) = fPos(xPos);
        NewHistCounts(xNeg+1,b,r) = flip(fNeg(xNeg));
    end
end
end

