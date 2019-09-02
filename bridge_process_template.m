%% extract and filter

%parameters
o = iss;
o.nRounds = 7;
o.nExtraRounds = 1;         %Treat Anchor channel as extra round
o.InputDirectory = 'C:\Users\...\Experiment1\raw_data';     %Folder path of raw data

%FileBase{r} is the file name of the raw data of round r in o.InputDirectory
o.FileBase = cell(o.nRounds+o.nExtraRounds,1);
o.FileBase{1} = 'round0';
o.FileBase{2} = 'round1';
...
o.FileBase{o.nRounds+o.nExtraRounds} = 'anchor';    %Make sure the last round is the anchor

o.RawFileExtension = '.nd2';
o.TileDirectory = 'C:\Users\...\Experiment1\tiles'; 
o.DapiChannel = 1;
o.ReferenceRound = o.nRounds+o.nExtraRounds;
o.OutputDirectory = 'C:\Users\...\Experiment1\output';  

%run code
o = o.extract_and_filter;
save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');

%% register

%parameters
o.TileSz = 2048;
o.AnchorChannel = 7;
o.RegMinSize = 8000; 
o.MaxRegShift = 50;
o.MaxOverlapFract = 0.2;

%run code
o = o.register;
save(fullfile(o.OutputDirectory, 'oRegister'), 'o', '-v7.3');

%% find spots

%parameters
o.SmoothSize = 2;
o.DetectionThresh = 300;
o.DetectionRadius=2;
o.IsolationRadius1 = 4;
o.IsolationRadius2 = 14;
o.nBP = 7;

%If a channel or round is faulty, you can ignore it by selecting only the
%good ones in o.UseChannels and o.UseRounds.
o.UseChannels = 1:o.nBP;
o.UseRounds = 1:o.nRounds;

o.FirstBaseChannel = 1;
o.InitialShiftChannel = 7;      %Channel to use to find initial shifts between rounds
o.PcDist = 3; 
o.MinPCMatches = 50; 
o.bpLabels = {'0', '1', '2', '3','4','5','6'}; %order of bases

%run code
o = o.find_spots;
save(fullfile(o.OutputDirectory, 'oFind_spots'), 'o', '-v7.3');

%% call spots

%parameters
%Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
%the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
o.CodeFile = 'C:\Users\...\Experiment1\codebook.txt';

%run code
o = o.call_spots;
save(fullfile(o.OutputDirectory, 'oCall_spots'), 'o', '-v7.3');

%% plot results

o.CombiQualThresh = 0.5;
BigDapiImage = imread(o.BigDapiFile);
figure(100);
o.plot(BigDapiImage);

%iss_view_codes(o,100,'[');

