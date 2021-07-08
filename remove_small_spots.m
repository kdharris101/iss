function Keep = remove_small_spots(Image,PeakPos)
%% Keep = remove_small_spots(Image,yPeak,xPeak)
%This takes peaks in found in Image at locations (yPeak,xPeak) then it
%exlcudes those peaks who have an adjacent pixel with negative intensity.
%Image: where peaks were found
%PeakPos: [yPeak, xPeak]

SzY = size(Image,1);
SzX = size(Image,2);
T = [[1,0];[0,1];[-1,0];[0,-1]];
Keep = zeros(size(PeakPos,1),size(T,1));
for i=1:size(T,1)
    ModPeakPos = PeakPos+T(i,:);
    ModPeakPos(ModPeakPos(:,1)>SzY,1) = SzY;
    ModPeakPos(ModPeakPos(:,2)>SzX,2) = SzX;
    ModPeakPos(ModPeakPos<1) = 1;
    Keep(:,i) = Image(sub2ind(size(Image),ModPeakPos(:,1),ModPeakPos(:,2)))>0;
end
Keep = all(Keep,2);
end

