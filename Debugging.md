# Debugging
This document gives some tips to overcome some common problems in the various steps of the in situ sequencing pipeline.

## [extract_and_filter.m](https://github.com/jduffield65/iss/blob/7x7-No-Anchor/%40iss/extract_and_filter.m)
**Problem 1**

The only real problem that can be encountered in this step is to do with the filtering. At the moment, we do a top hat filtering with the intention of extracting spots from the background. Below atre some examples of this step with vaying size of filters (filters are disks of radius r pixels).

<p float="left">
<img src="DebugImages/extract_and_filter/NoFilter.png" width = "350"> 
<em> No Filter</em>
<img src="DebugImages/extract_and_filter/r15.png" width = "350">
<em> r = 15</em>
</p>

<p float="left">
<img src="DebugImages/extract_and_filter/r6.png" width = "350"> 
<em> r = 6 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</em>
<img src="DebugImages/extract_and_filter/r2.png" width = "350">
<em> r = 2</em>
</p>

From these images, it can be seen that the top two images don't fully extract the spots. There are a few places where the spots can't be distinguished too well from the background. The last image, with r = 2, removes too much information and the spots look deformed. So a comprimise with r=6 seems to be the best.

**Solution**

[The filter used to do the top hat filtering is](https://github.com/jduffield65/iss/blob/0f30eecc859066c9a9d365229a2686aa1cfd2808/%40iss/extract_and_filter.m#L50-L58):
```matlab
SE = strel('disk', o.ExtractR);
```
As default, ```o.ExtractR = 'auto'``` which means ```o.ExtractR = 1/pixelsize```. This usually works but if you see spots in the raw data which are not emphasized by the filtering, there is the option to manually set ```o.ExtractR```. It should be approximately the radius of the spots in the raw data. 

**Problem 2**

Another potential problem stems from the fact that the datatype of the saved filtered tiles is ```uint16``` (integers in the range 0 to 65535). So if there are a lot of values in the filtered images above 65,535 or equally well below (rounding to integers, loses details in decimal places) then there will be a significant loss of information. You may notice this in the [spot detection step](https://github.com/jduffield65/iss/blob/0f30eecc859066c9a9d365229a2686aa1cfd2808/%40iss/find_spots.m#L54)  if there are a lot of spots very close to each other (within a small region, due to this rounding there may be lots of pixels with the same value and hence when finding the local maxima, multiple are found and identified as spots). 

**Solution**
Ideally, we want the multiply each image by the same amount so overall the pixel values fill the range 0-65,535. To fix this problem, [just before saving the tile](https://github.com/jduffield65/iss/blob/0f30eecc859066c9a9d365229a2686aa1cfd2808/%40iss/extract_and_filter.m#L98-L101) I would add the following:
```matlab
if t == 1 && r == 1 && c == 1
  ScaleFactor = 10000/max(max(IFS));
end
IFS = IFS*ScaleFactor;
```
The value of 10,000 is chosen as to give enough information without getting too close to the maximum (later images might have higher max values). I would say though, that I have never yet had this problem in 2D.
