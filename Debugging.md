# Debugging
This document gives some tips to overcome some common problems in the various steps of the in situ sequencing pipeline.

## [extract_and_filter.m](https://github.com/jduffield65/iss/blob/7x7-No-Anchor/%40iss/extract_and_filter.m)

For each image corresponding to a particular imaging round, colour channel and tile, this step first applies focus stacking. It then filters the resultant 2D image and saves the result as tiffs. 

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


## [register.m](https://github.com/jduffield65/iss/blob/7x7-No-Anchor/%40iss/register.m)

In the register step, the different tiles are stitched together. The aim is to [find all the shifts between neighbouring tiles](https://github.com/jduffield65/iss/blob/e054355308182990de1ef397377d6ab6b0e66768/%40iss/register.m#L58-L90), from which a global coordinate system can be created. Thus, the only problems in this step stem from finding the wrong shifts.

**Problem 1**

When you run ```o = o.register;```, all the shifts between tiles are read out. However, sometimes the shifts may be NaN:

<p float="left">
<img src="DebugImages/register/NaNShifts.png" width = "350"> 
</p>


**Solution**

This arises because the maximum correlation is below the threshold imposed by ```o.RegCorrThresh```. Too see the shifts found anyway, one should lower this: ```o.RegCorrThresh = -100;```, and in this case we get:

<p float="left">
<img src="DebugImages/register/FixedNaNShifts.png" width = "350"> 
</p>

However, the fact that the correlation of these shifts are below the threshold suggests that the shifts are not correct so the following problem is also likely to occur.

**Problem 2**

The shift values should follow approximately the pattern:

* right: [0,-1820]
* down: [-1820, 0]

So if a shift value is drastically different to this or if the cc value is significantly lower than the others (as with Tile 24 - right below), there is likely to be a problem.

<p float="left">
<img src="DebugImages/register/CorrelationList.png" width = "350"> 
</p>

To visualise the poroblem more clearly, you can set ```o.Graphics = 2;```. Then for each shift found, a plot of the Log Correlation will be displayed, with the best shift labelled by a cross. The plot for Tile 24 - right is below:

<p float="left">
<img src="DebugImages/register/LogCorrelation.png" width = "650"> 
</p>

Now we can clearly see the problem - the search over shifts was restricted in such a way that the local maxima was excluded.

**Solution**

The search is restricted as such:

* non-overlapping direction: ```-o.MaxRegShift < Shift < o.MaxRegShift```
* overlapping direction: ```-o.TileSz + o.MaxRegShift < Shift < -o.TileSz*(1-o.MaxOverlapFract)```

So changing the parameters of ```o.MaxRegShift``` and ```o.MaxOverlapFract``` are very significant for ```o.register```. The default values should of ```o.MaxRegShift = 50``` and ```o.MaxOverlapFract = 0.2``` should work well though. Indeed, if these values are used for the example shown above, we get:

```Tile 24 (3, 4), right: shift 6 -1816, cc 0.863064```

**Other tips**

* The first place to start is by looking at these correlation plots. They tell you if there is a peak where there should be that is just not being identified or if there is no peak at all. 
* If there is no peak at all, it may be worth playing around with the value of ```o.RegMinSize```. This controls the minimum size of overlap between tiles.


## [find_spots.m](https://github.com/jduffield65/iss/blob/7x7-No-Anchor/%40iss/find_spots.m)

This step first identifies all the spots in the anchor round. It then identifies all the spots in each round and colour channel. From these, point cloud registration is applied to find a transformation between the anchor round and every other round and colour channel for each tile. This transform is then applied to every spot in the anchor round to find its corresponding pixel value in each of the other rounds and colour channels. So for each spot, we end up with a vector of length ```o.nRounds*o.nBP```.

**Problem 1**

The first obstacle is identifying spots. A good way of visualsing spots [at this point](https://github.com/jduffield65/iss/blob/3757b0ac21a2b04e769ae101ce7d203bece3b809/%40iss/find_spots.m#L54) is by again setting ```o.Graphics = 2```. A plot like the one below would then appear for each tile in turn. 

<p float="left">
<img src="DebugImages/find_spots/AnchorSpots.png" width = "650"> 
</p>

The coloured circles highlight three potential problems:

* Green: may consider this to be a falsely identified spot.
* Cyan: may consider this to be a spot that was not identified.
* Yellow: may consider these spots to be too close together.

**Solution**

The spots are detected by [fnding the local maxima](https://github.com/jduffield65/iss/blob/3757b0ac21a2b04e769ae101ce7d203bece3b809/%40iss/detect_spots.m#L63) that are above the threshold set by ```o.DetectionThresh```.

The green and cyan problems can be addressed by changing ```o.DetectionThresh```. This has a default value of 300 (the value used in the image), however if this isn't adequate a good way of curating it is by looking at the detected spots plots and identifying the pixel values of spots falsely identified or neglected.

<p float="left">
<img src="DebugImages/find_spots/GreenSpot.png" height = "350"> 
<img src="DebugImages/find_spots/CyanSpot.png" height = "350">
</p>

From this we can see, that if we see the green spot as a problem, we should increase ```o.DetectionThresh``` above 365. Equally, if we see the cyan spot as an issue we should decrease ```o.DetectionThresh``` below 273.

The yellow problem can be addressed by changing ```o.DetectionRadius```. This controls the size of the neighbourhood in which the local maxima is found. So increasing this value increases the minimum allowed separation of spots. 

**Problem 2**

It is possible that when looking at the detected spots, you may see some spots that have pixel values below the threshold imposed by ```o.DetectionThresh```.

**Solution**

This situation arises because by default the threshold is [iteratively decreased until a certain number of spots are found](https://github.com/jduffield65/iss/blob/3757b0ac21a2b04e769ae101ce7d203bece3b809/%40iss/detect_spots.m#L55-L66). The minimum number of spots required is set by ```o.minPeaks``` and has a default value of 1000. For every iteration, where the number of spots found falls below this, the threshold decreases by ```o.ThreshParam``` which has a default value of 5. 

This should not really affect the anchor round as this is the round with the most spots. However, when [detecting spots in each round and colour channel](https://github.com/jduffield65/iss/blob/3757b0ac21a2b04e769ae101ce7d203bece3b809/%40iss/find_spots.m#L154), some images required a lower threshold than others. Hence this was introduced to get around this, and also around 1000 spots were required in each image for the [registration](https://github.com/jduffield65/iss/blob/3757b0ac21a2b04e769ae101ce7d203bece3b809/%40iss/find_spots.m#L170) to work. 

If this does cause problems though e.g. detecting lots of spots that aren't there, you can get rid of it by setting ```o.minPeaks = 1```.



**Other tips**

* It is easier to debug by just focussing on a particular tile. You can do this with ```o.EmptyTiles```. This is a nY x nX array, and all tiles with an index value of 0 are used. So to just use the tile with YX position (y,x) you can do the following:

```matlab
o.EmptyTiles = o.EmptyTiles+1;
o.EmptyTiles(y,x) = 0;
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Then if you do ```o = o.find_spots;```, it will do the whole section with just one tile.

* Similarly, you can select which rounds and colour channels to use using ```o.UseRounds``` and ```o.UseChannels``` respectively. Note in this instance, both the first and round and channel are indicated by 1. Usually the colour channels, in which the spots are clearest are 5,6,7 and to only use you write ```o.UseChannels = [5,6,7];```


