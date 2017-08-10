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
        InputDirectory;
        
        % where the top-hat image files are kept (one for each tile and round, all channels)
        TileDirectory;
        
        % where .mat file outputs and intermediates are kept
        OutputDirectory;
        
        % code file for spot calling
        CodeFile;
        
        %% parameters: stuff that might vary between experiments
        
        % one more round after the main rounds for Sst and Npy
        nExtraRounds = 1; 
        
        % For single-gene rounds, a n by 3 structure array of {GeneName, round, channel, threshold}
        ExtraCodes = {'Sst', 6, 3, 500; 'Npy', 6, 5, 400};
        
        % BasePair labels
        bpLabels = {'A', 'C', 'G', 'T'};
        
        %% parameters: registration and alignment
        
        % correlation threshold for image alignment. Can be low since 
        % chance correls are miniscule. But shouldn't be too low since 
        % microscopes can have uneven pixel intensities, generating
        % spurious correlations of 0 offset
        RegCorrThresh = [.2 .6]; 
                
        % minimum size overlap for tile matching (pixels - so side squared)
        RegMinSize = 100^2; 
        
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
        
        %% parameters: spot calling
        % normalizes spot fluorescence so this percentile = 1
        SpotNormPrctile = 98;
        
        % score and intensity thresholds to plot a spot (combi codes)
        CombiQualThresh = .8; 
        CombiIntensityThresh = .1;
        
        %% parameters: cell segmentation
        
        % percentile threshold for dapi image (after imadjust)
        DapiThresh = 90;
        
        % how close local maxima can be together (in pixels)
        DapiMinSep = 7;
        
        % how big cell radius needs to be to be detected (in pixels)
        DapiMinSize = 5;
        
        % how much margin around each cell to assign to it (in pixels)
        DapiMargin = 10;
        
        % minimum area in pixels for a cell not to be dropped
        MinCellArea = 200;
        
        %% parameters: cell calling
        
        % CellCallRegionXY(:,1:2): XY coords of polygod outlining region to
        % look for cells in
        CellCallRegionYX;
       
        % how many neighboring cells to consider for each spot
        nNeighbors = 3; 
        
        % gamma dist prior for each spot
        rSpot = 2; 
        
        % gamma dist shape for scaling whole gene
        rGene = 20; 
        
        % probability of misread per pixel
        MisreadDensity = 1e-5; 
        
        % prior on how much to scale scRNAseq down by to match iss
        Inefficiency = .2; 
        
        % how much to add to expression so as not to get infities
        SpotReg = .1; 
        
        % additional log likelihood for a gene inside the segmented region
        InsideCellBonus = 2; 
        
        % converges when no probabilities have changed more than this
        CellCallTolerance = .02; 
        
        % converges when no probabilities have changed more than this
        CellCallMaxIter = 100; 

        
        %% parameters: stuff that should be the same between experiments
        
        % for dapi images: scale is downsampling factor for final .fig file
        DapiChannel = 1;
        
        % which channel of each file is anchor images
        AnchorChannel = 2;
        
        % which sequencing round to align all others to
        ReferenceRound = 2;
        
        % how many combinatorial sequencing rounds 
        nRounds = 5;
        
        % Number of possible basepairs (always 4 for life as we know it but
        % hardcoding is bad programming style!)
        nBP = 4;
        
        % tile size (assumed square tiles)
        TileSz = 2048;
        
        % expected step size (in coordinates returned by bioformats)
        MicroscopeStepSize = 2048;
        
        %% variables: filenames
        
        % TileFiles{r, y, x} is full pathname of tile (y,x) on round r
        % (mutlicolor tiff file). Empty if no file there
        TileFiles; 
        
        % EmptyTiles{y,x} is zero if that tile location wasn't scanned
        EmptyTiles;
        
        % FileBase{r} is local filename of round r .czi file; other
        % extensions added in output files
        FileBase;
        
        % BigDapiFile is full path of stitched DAPI image (tiff format)
        BigDapiFile;
        
        % CellMapFile is full path of cell map (.mat format)
        CellMapFile;
        
        %% variables: registration
        
        % TilePoxYX(t): grid position of tile t in integers. t is what you
        % find in the file name - so it only counts non-empty tiles
        TilePosYX;
        
        % RefPos(y, x): origin of tile (y,x) in pixels on
        % reference round relative to global coordinate frame
        RefPos;
        
        % RelativePos(r, 1:2, t1, t2): origin of tile t2 on 
        % reference round minus origin of tile t1 round r. In other words,
        % Im2(x;rr) = Im1(x + RelativePos; r). Nan if not defined
        % t1 and t2 are linear indices for the tile (y,x)
        RelativePos; 
        
        % RegistrationCorr(r,t1,t2): image correlation for the anchor channel
        % on each registration
        RegistrationCorr;
        %% variables: spot calling outputs
       
        % cSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        cSpotColors;
        
        % cSpotIsolated(Spot) is a binary array saying if the spot is well isolated
        % again for combinatorial spots only
        cSpotIsolated; 
        
        % SpotCombi: binary array which is one for combinatorial spots,
        % zero for extra spots
        SpotCombi;
        
        % SpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        SpotGlobalYX;
        
        % SpotCodeNo(Spot): code number for each spot (uint16). Both combinatorial and extra spots
        SpotCodeNo;
        
     
        % SpotScore(Spot): score saying how well the code fits (0...1).
        % 1 for extras
        SpotScore;
        
        % SpotIntensity(Spot): RMS intensity of the spot. Zero for
        % extras!
        SpotIntensity;
        
        % CharCodes{Code}: text representation of each code in the
        % codebook. For extra spots, this just says "EXTRA"
        CharCodes;
        
        % GeneNames{Code}: gene name for each code. Note that some genes
        % can have multiple entries, because they have multiple codes or
        % because they have extra rounds
        GeneNames;

		
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

 





