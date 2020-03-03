# How to run
The only file that you need to run to obtain and save the data is [bridge_process_template.m](https://github.com/jduffield65/iss/blob/master/bridge_process_template.m). The following will explain the changes to this file that need to made in order for it work with your data.

## File names
There are another of file/folder paths which need to be given:
* [```o.InputDirectory```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L7): This is the path of the folder that contains the raw data, of type specified by [```o.RawFileExtension```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L20), e.g. ```.nd2```. An example of what this folder typically looks like is given below:

<p float="left">
<img src="DebugImages/README/InputDirectory.png" width = "450"> 
</p>

* [```o.FileBase```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L10-L18): These are the names of the files within ```o.InputDirectory``` (minus the extension). For the above example, we would set
<pre>
o.FileBase{1} = 'Exp1_r0';
o.FileBase{2} = 'Exp1_r1';
&#8942
o.FileBase{7} = 'Exp1_r6';
o.FileBase{8} = 'Exp1_anchor';
</pre>
You need to make sure that ```o.FileBase{```[```o.ReferenceRound```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L24)```}``` is set to the anchor round. Also, the other rounds must be in the correct imaging order.

* [```o.TileDirectory```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L21): This is the path for the folder that you would like the filtered images for each tile, round and colour channel to be saved to. The file named as  ```o.FileBase{r}_tT.tif``` contains all the colour channel images for round r, tile T.

* [```o.OutputDirectory```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L26): This is the path for the folder that you would like the iss objects after each step of the pipeline to be saved.

* [```o.CodeFile```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L105): This is the path for the file containing the code for each gene. The file should be a text file containing two columns, the first being the gene name. The second is the code of length [```o.nRounds```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L5), containing numbers in the range from 0 to [```o.nBP```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L74)-1. An example codebook for ```o.nRounds=o.nBP=7``` is shown below:

<p float="left">
<img src="DebugImages/README/CodeBook.png" width = "250"> 
</p>

## Stitching and registration parameters
These are parameters that slightly affect how the stitching of tiles and registration between rounds and colour channels works. The default values in [bridge_process_template.m](https://github.com/jduffield65/iss/blob/master/bridge_process_template.m) should work most of the time but there are some cases when they may need altering.

* [```o.InitialShiftChannel```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L82): This is the channel used to register all rounds to the anchor, so ensure it is one of the best colour channels. 

* [```o.RegSearch```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L61-L64): The algorithm for stitching together tiles only looks at shifts in the range specified by ```o.RegSearch```. The values that may need changing here are ```o.RegSearch.South.Y``` and ```o.RegSearch.East.X```. The default values of ```-1900:o.RegStep(2):-1700``` are heavily dependent on [```o.TileSz```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L44). I.e. we only consider overlaps between 148 - 348 pixels which covered the expected value of about 10% (205 pixels). If the expected overlap or [```o.TileSz```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L44) is different to this though, these values will need changing. 

* [```o.FindSpotsSearch```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L88-L89): The algorithm for finding shifts between the anchor round and each imaging round only looks at shifts in the range specified by ```o.FindSpotsSearch```. The default values assume a relatively small shift between of absolute value less than 100 pixels. However, if you know one round had a particularly large shift, you will need to change this range. You can also specify a different range for each round through:
```matlab
o.FindSpotsSearch = cell(o.nRounds,1);
for r = o.UseRounds
    o.FindSpotsSearch{r}.Y = MinYShift:o.FindSpotsStep(1):MaxYShift;
    o.FindSpotsSearch{r}.X = MinXShift:o.FindSpotsStep(1):MaxXShift;
end
```

## Visualising results
To visualise the results, load in the final saved iss object which should be named [```oCallSpots.mat```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L111). Then, load in the  background dapi image and run [o.plot](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L116-L117). This will show you the gene assignments given by the file [```o.call_spots```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L109) which is achieved by taking the dot product of the normalised spot and bleed codes. Only the results where this dot product is above [```o.CombiQualThresh```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L115) are shown. An example plot is given below with ```o.CombiQualThresh = 0.7```.

To change the value of the threshold, simply set ```o.CombiQualThresh = NewValue``` and then run ```iss_change_plot(o)```. The above data with ```o.CombiQualThresh = 0.9``` is shown below:

With [```o.CallSpotsCodeNorm = 'WholeCode'```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L108), each spot and gene code has L2 norm of 1 so the maximum value of the dot product hence ```o.CombiQualThresh``` is 1 (Recommend ```o.CombiQualThresh ~ 0.7```). However, with ```o.CallSpotsCodeNorm = 'Round'```, each round in each code has L2 norm of 1 so each code has L2 norm of [```o.nRounds```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L5). So in this case, the max value of the dot product hence ```o.CombiQualThresh``` would be ```o.nRounds``` (Recommend ```o.CombiQualThresh ~ 4``` for ```o.nRounds=7```. Note that if you change ```o.CallSpotsCodeNorm```, then you need to run [```o.call_spots```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L109) again. The justification for setting ```o.CallSpotsCodeNorm = 'Round'``` is that with no bleed through, we expect each spot to appear in one colour channel in each round so we want to give each round equal weighting, but with ```o.CallSpotsCodeNorm = 'WholeCode'```, a particularly intense round would dominate the others.
