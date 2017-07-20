function [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] =...
    get_ome_tilepos(reader)
% [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] =...
%     get_ome_tilepos(reader)
% input: OMERO reader or microscopy image file that contains OMERO metadata
% Xiaoyan, 2017


% get reader if input is file name (character)
if ischar(reader)
    reader = bfGetReader(reader);
end

nSeries = reader.getSeriesCount();
nChannels = reader.getSizeC;

omeMeta = reader.getMetadataStore();
nZstacks = omeMeta.getPixelsSizeZ(0).getValue();
pixelsize = double(omeMeta.getPixelsPhysicalSizeX(0).value());

xypos  = [];
for i = 1:nSeries
    try
        posX = double(omeMeta.getPlanePositionX(i-1, 0).value());
        posY = double(omeMeta.getPlanePositionY(i-1, 0).value());
        xypos = [xypos; posX, posY];
    end
end

xypos = bsxfun(@minus, xypos, min(xypos,[],1))/pixelsize;
nSerieswPos = size(xypos,1);

% close
reader.close();

end
