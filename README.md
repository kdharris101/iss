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
<img src="DebugImages/README/CodeBook.png" width = "450"> 
</p>
