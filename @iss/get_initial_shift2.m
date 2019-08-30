function [BestShift,BestScore] = get_initial_shift2(o, y, x0,Search,section)
% Finds the initial shifts to give to PCR algorithm between rounds for each
% tile. Does this by finding the colour channel with the most spots for
% each round. Then uses ImRegFft2 to find the shift between this base binary
% image and the anchor channel binary image
%
% inputs:
% y is a cell containing the centered YX location of all spots in round r 
% , colour channel c for all tiles in units of XY pixels in find_spots. In
% register, y is equivalent to x0 but for another tile.
%
% x0 is a cell containing the non centered YX location of spots in the 
% anchor channel for the current tile.
%
% Search.Y,Search.X are the ranges in XY of shifts to search
%
% section specifies which part of the pipeline we are on: Register or
% FindSpots
%

%o = o.get_initial_shift(AllBaseLocalYXZ{t,b,r},RawLocalYXZ{t})
%D0 = zeros(nTiles,3,o.nRounds);

%y = AllBaseLocalYXZ{t,b,r};
%x0 = RawLocalYXZ{t};

if strcmpi(section, 'FindSpots')
    %centre anchor channel spots
    x = x0 - [o.TileSz/2,o.TileSz/2];
elseif strcmpi(section, 'Register')
    x = x0;
else
    warning('Have not specified which part of pipeline is being run')
end

% make kd tree - default options!
k = KDTreeSearcher(y);

%Find all permutations of shifts to try
%TO SPEED UP MAYBE DO AN INITIAL WIDE SEARCH WITH WIDER SPACING in X,Y???

[A,B] = meshgrid(Search.Y,Search.X);
c=cat(2,A,B);
shifts=reshape(c,[],2);

%Error = zeros(size(shifts,1),1);
Score = zeros(size(shifts,1),1);
%parfor i = 1:size(shifts,1)
for i = 1:size(shifts,1)
    xShifted = x+shifts(i,:);
    [~,Dist] = k.knnsearch(xShifted);
    %UseMe = Dist<1000;   
    %Error(i) = mean(Dist(UseMe>0).^2);
    %Error(i) = mean(Dist.^2);
    Score(i) = sum(exp(-Dist.^2/(2*o.ShiftScoreThresh^2)));
end
%min(Error)
BestScore = max(Score);
BestShift = shifts(Score == BestScore,:);

if o.Graphics == 2
    shiftsNew = shifts-min(shifts)+1;
    
    ScoreImage = zeros(max(shiftsNew));
    for i = 1:size(shiftsNew,1)
        ScoreImage(shiftsNew(i,1),shiftsNew(i,2)) = Score(i);
    end
    
    figure(83452);
    xIndex = min(shifts(:,2)):max(shifts(:,2));
    yIndex = min(shifts(:,1)):max(shifts(:,1));
    imagesc('XData',xIndex,'YData',yIndex,'CData',ScoreImage);
    
    xlim([min(shifts(:,2)) max(shifts(:,2))]);
    ylim([min(shifts(:,1)) max(shifts(:,1))]);
    xlabel('X');
    ylabel('Y');
    
    hold on
    scatter(BestShift(2),BestShift(1),80,'kx');
    hold off
end


fprintf('\n');
