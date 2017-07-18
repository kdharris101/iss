function o = segment_dapi(o, Dapi)
% CellMap = iss_segment_dapi(DapiIm)
%
% segments a DAPI image and assigns each pixel to a cell. 
%
% Output CellMap is same size as input DapiIm, with integer entries for
% each pixel, saying which cell it belongs to. (Zero if unassigned)

Dapi = imadjust(Dapi); % contrast enhancement
ImSz = size(Dapi);
Debug = 0;
%% threshold the map
ThreshVal = prctile(Dapi(:), o.DapiThresh);

bwDapi = imerode(imfill(Dapi>ThreshVal, 'holes'), strel('disk', 2));

if Debug
    figure(300)
    subplot(2,1,1);
    imagesc(Dapi); 
    subplot(2,1,2);
    imagesc(bwDapi);
    colormap bone
end
%% find local maxima 
dist = bwdist(~bwDapi);
dist0 = dist;
dist0(dist<o.DapiMinSize)=0;
ddist = imdilate(dist0, strel('disk', o.DapiMinSep));
clear dist dist0
impim = imimposemin(-dist0, imregionalmax(ddist));

if Debug
    figure(301);
    subplot(2,2,1)

    imagesc(dist);
    subplot(2,2,2)
    subplot(2,2,3)
    imagesc(maxxes);
    subplot(2,2,4)
    imagesc(impim);
end
%% segment
% remove pixels at watershed boundaries
bwDapi0 = bwDapi;
bwDapi0(watershed(impim)==0)=0;

% assign all pixels a label
labels = uint32(bwlabel(bwDapi0));
[d, idx] = bwdist(bwDapi0);

% now expand the regions by a margin
CellMap = zeros(ImSz, 'uint32');
Expansions = (d<o.DapiMargin);
CellMap(Expansions) = labels(idx(Expansions));

if Debug
    figure(302)
    subplot(2,2,1);
    image(label2rgb(ws, 'jet', 'w', 'shuffle'));

    subplot(2,2,2);
    imagesc(bwDapiSep);

    subplot(2,2,3);
    imagesc(labels);

    % now give every pixel to its nearest neighbor
    subplot(2,2,4);
    image(colors);
end
%% make image with boundaries
Boundaries = (CellMap ~= imdilate(CellMap,strel('disk', 1)));
DapiBoundaries = Dapi;
DapiBoundaries(Boundaries) = .3 * max(Dapi(:));

imwrite(DapiBoundaries, [o.OutputDirectory '\background_boundaries.tif']);
save CellMap CellMap
% figure(904375)
% imagesc(DapiBoundaries)


return

%% make a pseudocolored image

colors = label2rgb(CellMap, 'hsv', 'w', 'shuffle');

NormDapi = min(double(Dapi)/double(max(Dapi(:))),1);
pColorIm = bsxfun(@times, double(colors)/255, NormDapi);
Boundaries = (CellMap ~= imdilate(CellMap,strel('disk', 1)));
[y, x] = find(Boundaries);
for i=1:3
    pColorIm(sub2ind(size(pColorIm), y, x, i*ones(size(y))))=.3; %ones(size(y))*i]))=1;
end

figure(904375)
image(pColorIm)

%imdilate(bwDapiSep, strel('disk', 1))-bwDapiSep;
%pColorIm

return
% subplot(2,2,3)
% Props = regionprops(bwDapiSep, 'EquivDiameter');
% histogram([Props.EquivDiameter],0:max([Props.EquivDiameter]));
% 
% 
% excludeDapi = cellfun(@(v) sum(maximaDapi(v))==0, propDapi(2,:));       % no maxima
% excludeDapi = excludeDapi | ~cellfun(@(v) v>8 & v<100, propDapi(1,:));  % size threshold
%%
%Thresh = 500;%prctile(Dapi(:), 90);
Thresh = graythresh(Dapi)*max(Dapi(:))
bim = (Dapi>Thresh);
imagesc(bim);
colormap bone
%%
figure(302)
bwdim = bwdist(~bim);
imagesc(bwdim);
colorbar
%%
% figure(303)
% Outside = find(Dapi<Thresh);
% gdim = graydist(Dapi, Outside);
% imagesc(gdim)
%%
figure(304)
imax = imregionalmax(bwdim);
imagesc(imax.*bwdim)
%%
figure(305)
bwdim2 = bwdim; 
bwdim2(imax) = bwdim2(imax)-1;
imax2 = imregionalmax(bwdim2);
imagesc(imax2.*bwdim2);
%% 
figure(306);
iimp = imimposemin(-bwdim, imax2 & (bwdim>9));
imagesc(iimp);
%%
L = watershed(iimp);
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
image(bsxfun(@times, double(Lrgb), double(bim)));

%%
se = strel('disk', 4);
%grad = imdilate(Dapi, se) - imerode(Dapi, se);
hy = fspecial('sobel');
gy = double(imfilter(Dapi, hy, 'replicate'));
gx = double(imfilter(Dapi, hy', 'replicate'));
grad = sqrt(gy.^2 + gx.^2);
imagesc(grad)