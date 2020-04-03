function o = get_extract_scale(o,nChannels,nZstacks,bfreader,SE,DapiSE,r,t)
%This finds the scale by which to multiply each filtered image. This is
%based on the max intensity of tile t round r across all color channels.
%It also writes the images for this tile and round to o.TileDirectory

if r~=o.AnchorRound; fprintf('Obtaining ExtractScale for imaging rounds... ');
else; fprintf('Obtaining ExtractScale for anchor round... ');end
                
MaxPixelValue = 0;
RawImage_all = zeros(o.TileSz,o.TileSz,nChannels);
IFS_all = zeros(o.TileSz,o.TileSz,nChannels);
for c=1:nChannels
    I = cell(nZstacks,1);
    for z = 1:nZstacks
        iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
        I{z} = bfGetPlane(bfreader, iPlane);
    end
    
    % focus stacking
    I_mod = o.fstack_modified(I);
    RawImage_all(:,:,c) = I_mod;
    if c == o.DapiChannel && r == o.AnchorRound
        IFS_all(:,:,c) = imtophat(I_mod, DapiSE);
    else
        I_mod = single(padarray(I_mod,(size(SE)-1)/2,'replicate','both'));
        IFS = convn(I_mod,SE,'valid');
        if max(IFS(:))>MaxPixelValue
            MaxPixelValue = max(IFS(:));
        end
        IFS_all(:,:,c) = IFS;
        clearvars I_mod I IFS %Free up memory
    end
end

if r==o.AnchorRound
    o.ExtractScaleAnchor = o.ExtractScaleNorm/MaxPixelValue;
    IFS_all(:,:,setdiff(1:nChannels,o.DapiChannel)) = ...
        IFS_all(:,:,setdiff(1:nChannels,o.DapiChannel))*o.ExtractScaleAnchor;
else
    o.ExtractScale = o.ExtractScaleNorm/MaxPixelValue;
    IFS_all = IFS_all*o.ExtractScale;
end

if r~=o.AnchorRound; fprintf('Value is %.2f\n', o.ExtractScale);
else; fprintf('Value is %.2f\n', o.ExtractScaleAnchor);end

if o.Graphics
    o.TilePlot(IFS_all,RawImage_all,r,nChannels);
end

for c=1:nChannels
    if c == o.DapiChannel && r == o.AnchorRound
        IFS = uint16(IFS_all(:,:,c));
    else        
        if c ~= o.AnchorChannel && r == o.AnchorRound
            IFS = uint16(IFS_all(:,:,c)+o.TilePixelValueShift);
        else
            %Determine auto thresholds
            IFS = IFS_all(:,:,c);
            o.AutoThresh(t,c,r) = median(abs(IFS(:)))*o.AutoThreshMultiplier;
            if ismember(r,1:o.nRounds)
                %Get histogram data
                IFS = int32(IFS);
                o.HistCounts(:,c,r) = o.HistCounts(:,c,r)+histc(IFS(:),o.HistValues);
            end
            IFS = uint16(IFS+o.TilePixelValueShift);
        end
        
    end
    imwrite(IFS,...
        fullfile(o.TileDirectory,...
        [o.FileBase{r}, '_t', num2str(t), '.tif']),...
        'tiff', 'writemode', 'append');
    fprintf('Round %d tile %d colour channel %d finished.\n', r, t, c); 
end


end
