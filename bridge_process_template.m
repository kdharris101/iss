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
o.FileBase{3} = 'round2';
o.FileBase{4} = 'round3';
o.FileBase{5} = 'round4';
o.FileBase{6} = 'round5';
o.FileBase{7} = 'round6';
o.FileBase{8} = 'anchor';    %Make sure the last round is the anchor

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

%Anchor spots are detected in register2
o.DetectionRadius = 2;
o.SmoothSize = 2;     
o.IsolationRadius1 = 4;
o.IsolationRadius2 = 14;

o.DetectionThresh = 'medianx10';       %SEEMS TO BE THE BEST AUTO METHOD AT THE MOMENT
o.MinThresh = 20;
o.minPeaks = 1000;

%paramaters to find shifts between overlapping tiles
o.RegMinScore = 30;     
o.RegStep = [5,5];
o.RegSearch.South.Y = -1900:o.RegStep(1):-1700;
o.RegSearch.South.X = -50:o.RegStep(2):50;
o.RegSearch.East.Y = -50:o.RegStep(1):50;
o.RegSearch.East.X = -1900:o.RegStep(2):-1700;
o.RegWidenSearch = [50,50]; 

%run code
o = o.register2;
save(fullfile(o.OutputDirectory, 'oRegister'), 'o', '-v7.3');

%% find spots

%parameters
o.nBP = 7;

%If a channel or round is faulty, you can ignore it by selecting only the
%good ones in o.UseChannels and o.UseRounds.
o.UseChannels = 1:o.nBP;
o.UseRounds = 1:o.nRounds;

o.FirstBaseChannel = 1;

%Search paramaters
o.InitialShiftChannel = 7;      %Channel to use to find initial shifts between rounds
o.FindSpotsMinScore = 70;
o.FindSpotsStep = [5,5];
%FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
%structs - have a different range for each round: 
%o.FindSpotsSearch = cell(o.nRounds,1);
o.FindSpotsSearch.Y = -100:o.FindSpotsStep(1):100;
o.FindSpotsSearch.X = -100:o.FindSpotsStep(2):100;
%Make WidenSearch larger if you think you have a large shift between rounds
o.FindSpotsWidenSearch = [50,50]; 

o.PcDist = 3; 
o.MinPCMatches = 1;    %HACK SO IT GETS TO THE END
o.bpLabels = {'0', '1', '2', '3','4','5','6'}; %order of bases

%run code
o = o.find_spots2;
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

%iss_view_codes(o,100);

