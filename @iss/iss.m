%% iss: code for processing of in situ sequencing
% Kenneth D. Harris and Xiaoyan Qian
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% o = iss; % create structure, default parameters
% % change any parameters you want, and set o.FileN
% o = o.extract_and_filter; % create top-hat filtered tiffs for each tile
% o = o.find_spots; % find spot positions in global coordinates
% o = o.call_spots; % allocate each spot to a gene
% o = o.call_cells; % identify cells
%
% this current file iss.m just contains default parameter values

classdef iss
    properties
        %% interactive graphics mode. 1 means some, 2 means a lot.
        Graphics = 1;
        
        % region to show during cell calling
        CellCallShowCenter = [385 1100];
        CellCallShowRad = 200;
        ExampleCellCenter = [1670 250];

        
        %% parameters: file locations
        % where the input czi files are kept
        InputDirectory = '\\basket.cortexlab.net\data\kenneth\iss\170315_161220KI_4-1\Raw';
        
        % where the top-hat image files are kept (one for each tile and round, all channels)
        TileDirectory = '\\basket.cortexlab.net\data\kenneth\iss\170315_161220KI_4-1\FS_tophat_stack';
        
        % where .mat file outputs and intermediates are kept
        OutputDirectory = '\\basket.cortexlab.net\data\kenneth\iss\170315_161220KI_4-1\output';
        
        % code file for spot calling
        CodeFile = ['A:\Dropbox\Dropbox (Neuropixels)\170315_161220KI_4-1' '\codebook_unique.csv'];
        
        %% parameters: stuff that might vary between experiments
        
        % one more round after the main rounds for Sst and Npy
        nExtraRounds = 1; 
        
        % For single-gene rounds, a n by 3 structure array of {GeneName, round, channel}
        ExtraCodes = {'Npy', 6, 3; 'Sst', 6, 4};
        
        % BasePair labels
        bpLabels = 'TGCA';
        
        %% parameters: registration and alignment
        
        % correlation threshold for image alignment. Can be low since chance correls are miniscule
        CorrThresh = .2; 
        
        % minimum size overlap for tile matching
        MinSize = 100; 
        
        % distance scale for point cloud registration (pixels)
        PcDist = 3; 
        
        
        %% parameters: spot detection
        
        % smooth images before detecting fluorescence with a disk of this radius:
        SmoothSize = 1;
        
        % to detect spot, pixel needs to be above dilation with this radius
        DetectionRadius = 1;
        
        % and pixel raw fluorescence needs to be above this value:
        DetectionThresh = 300;
        
        % find isolated spots by annular filtering with these radii
        IsolationRadius1 = 2;
        IsolationRadius2 = 7;
        
        % annular filtered value needs to be less than this:
        IsolationThresh = 60;
        
        % for visualization during spot detection
        FindSpotsRoi = [1742 1755 213 227];

        %% parameters: cell segmentation
        
        % percentile threshold for dapi image (after imadjust)
        DapiThresh = 80;
        
        % how close local maxima can be together (in pixels)
        DapiMinSep = 7;
        
        % how big cell radius needs to be to be detected (in pixels)
        DapiMinSize = 5;
        
        % how much margin around each cell to assign to it (in pixels)
        DapiMargin = 10;
        
        % minimum area in pixels for a cell not to be dropped
        MinCellArea = 200;
        
        %% parameters: cell calling
       
        % how many neighboring cells to consider for each spot
        nNeighbors = 3; 
        
        % gamma dist prior for each spot
        rSpot = 2; 
        
        % gamma dist shape for scaling whole gene
        rGene = 100; 
        
        % probability of misread per pixel
        MisreadDensity = 1e-4; 
        
        % prior on how much to scale scRNAseq down by to match iss
        Inefficiency = .2; 
        
        % how much to add to expression so as not to get infities
        SpotReg = .2; 
        
        % additional log likelihood for a gene inside the segmented region
        InsideCellBonus = 3; 
        
        % converges when no probabilities have changed more than this
        CellCallTolerance = .01; 
        
        % converges when no probabilities have changed more than this
        CellCallMaxIter = 100; 

        
        %% parameters: stuff that should be the same between experiments
        
        % for dapi images: scale is downsampling factor for final .fig file
        DapiChannel = 1;
        
        % which channel of each file is anchor images
        AnchorChannel = 2;
        
        % which sequencing round to align all others to
        ReferenceRound = 2;
        
        % how many sequencing rounds (not counting single-gene rounds)
        nRounds = 5;
        
        % Number of possible basepairs (always 4 for life on earth)
        nBP = 4;
        
        % tile size (assumed square tiles)
        TileSz = 2048;
        
        % expected step size (in coordinates returned by bioformats)
        MicroscopeStepSize = 2048;
        
        %% variables: filenames
        
        % TileFiles{r, y, x} is full pathname of tile (y,x) on round r
        % (mutlicolor tiff file). Empty if no file there
        TileFiles; 
        
        % FileBase{r} is local filename of round r .czi file; other
        % extensions added in output files
        FileBase;
        
        %% variables: spot calling outputs
        % SpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system
        SpotGlobalYX;
        
        % SpotColors(Spot, Base, Round) contains spot color on each base
        % and round
        SpotColors;
        
        % SpotIsolated(Spot) is a binary array saying if the spot is well isolated
        SpotIsolated; 
        
        % SpotGene(Spot): number of gene for each spot
        SpotGene;
        
        % GeneName{Gene}: name of each gene
        GeneName;
        
        % SpotCode{Spot}: text representation of code for each spot 
        SpotCode;
        
        % SpotScore(Spot): score saying how well the code fits (0...1)
        SpotScore;
        
        % SpotIntensity(Spot): RMS intensity of the spot
        SpotIntensity;
        
        %% variables: cell calling outputs
        % pCellClass(cell, class); % prob each cell goes to each class: last class is zero expression
        pCellClass;
        
        % ClassNames(class): name of each class, taken from gSet file
        ClassNames;
        
        % pSpotCell(spot, cell): sparse array containing prob of each spot
        % to belong to each cell. Last cell is background
        pSpotCell;
        
    end
    
end

 





