# How to run
The only file that you need to run to obtain and save the data is [bridge_process_template.m](https://github.com/jduffield65/iss/blob/master/bridge_process_template.m). The following will explain the changes to this file that need to made in order for it work with your data.

## File names
There are another of file/folder paths which need to be given:
* [```o.InputDirectory```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L7): This is the path of the folder that contains the raw data, of type specified by [```o.RawFileExtension```](https://github.com/jduffield65/iss/blob/eb6d7c23acf2b59a18903511b25b34ecd756c05b/bridge_process_template.m#L20), e.g. ```.nd2```. An example of what this folder typically looks like is given below:

<p float="left">
<img src="DebugImages/README/InputDirectory.png" width = "450"> 
</p>




