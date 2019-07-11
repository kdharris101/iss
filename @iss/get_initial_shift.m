function o = get_initial_shift(o, y, x0, nTiles, MinSize)
% Finds the initial shifts to give to PCR algorithm between rounds for each
% tile. Does this by finding the colour channel with the most spots for
% each round. Then uses ImRegFft2 to find the shift between this base binary
% image and the anchor channel binary image
%
% inputs:
% y is a cell containig the centered YX location of all spots in all rounds 
% and colour channels for all tiles
%
% x0 is a cell containing the non centered YX location of spots in the 
% anchor channel for all tiles
%
% MinSize is for ImRegFft2. Needs to be much larger for these binary images
% than for the registration stage. Not sure by how much though
%

o.D0 = zeros(nTiles,2,o.nRounds);

fprintf('\nGetting initial shifts for tile   ');
for t=1:nTiles
    if o.EmptyTiles(t); continue; end    
    if t<10
        fprintf('\b%d', t);
    else
        fprintf('\b\b%d', t);
    end    
    
    %Make binary anchor channel image
    LocalImBinary = zeros(o.TileSz,o.TileSz);
    idx = sub2ind(size(LocalImBinary),x0{t}(:,1),x0{t}(:,2));
    LocalImBinary(idx) = 1;
    
    %Find color channel with most spots for each tile,round
    for r=o.UseRounds
        MostSpots = 0;
        BestChannel = 0;
        for b=o.UseChannels  
            if size(y{t,b,r},1) > MostSpots
                MostSpots = size(y{t,b,r},1);
                BestChannel = b;
            end
        end
        
        %Make binary image for this color channel
        BaseImBinary = zeros(o.TileSz,o.TileSz);
        BaseYX = y{t,BestChannel,r} + [o.TileSz/2,o.TileSz/2];
        idx2 = sub2ind(size(BaseImBinary),BaseYX(:,1),BaseYX(:,2));
        BaseImBinary(idx2) = 1;
        
        [o.D0(t,:,r), ~] = ImRegFft2(BaseImBinary,LocalImBinary, 0, MinSize);     
    end
end
fprintf('\n');
