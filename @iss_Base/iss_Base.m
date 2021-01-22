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

classdef iss_Base
    properties
        %%LogToFile = 1 if you want command window output to file
        LogToFile = 1;
        
        % interactive graphics mode. 1 means some, 2 means a lot.
        Graphics = 1;
        
        % region to show during cell calling
        CellCallShowCenter = [385 1100];
        CellCallShowRad = 200;
        
        % example cell to diagnose calling, plus class pair to compare (two
        % strings)
        ExampleCellCenter = [1670 250];
        CellCallDiagnosisPair = [];
        
        %% parameters: file locations
        %LogFile is where command window log output to
        LogFile;
        
        % where the input czi files are kept
        InputDirectory;
        
        % where the top-hat image files are kept (one for each tile and round, all channels)
        TileDirectory;
        
        % where .mat file outputs and intermediates are kept
        OutputDirectory;
        
        % code file for spot calling
        CodeFile;
        
        %% parameters: stuff that might vary between experiments
        
        % one more round after the main rounds for Anchor
        nExtraRounds = 1; 
        
        % For single-gene rounds, a n by 3 structure array of {GeneName, round, channel, threshold}
        ExtraCodes = {'Sst', 6, 3, 500; 'Npy', 6, 5, 400};
        
        % BasePair labels
        bpLabels = {'A', 'C', 'G', 'T'};
        
        %% parameters: extract_and_filter
        % To extract spots, a filter is used with inner radius ExtractR1 within 
        % which it is positive and outer radius ExtractR2 so the annulus between 
        % R1 and R2 is negative. Overall sums up to 0. R1 should be the
        % approximate radius of spot. R2 should be such that between 
        % R1 and R2 there is a dark region.
        ExtractR1 = 'auto';
        ExtractR2 = 'auto';
        
        % Each filtered image is multiplied by ExtractScale. This is
        % because the image is saved as uint16 so to gain information from
        % the decimal points, should multiply image so max pixel number is
        % in the 10,000s (less than 65,536). If 'auto', it sets to
        % ExtractScaleNorm/max(Tile ExtractScaleTile, round 1).
        % If ExtractScaleTile is 'auto', the default will be the middle tile.
        ExtractScaleNorm = 40000;   %Max pixel value in saved tile 1 round 1 will be
        %ExtractScaleNorm+TilePixelValueShift. Make sure this is about 10,000 below 65,536
        ExtractScaleTile = 'auto';
        ExtractScale = 'auto';
        ExtractScaleAnchor = 'auto';        %Computed in same way, but anchor can have different scaling
        
        % TilePixelValueShift is added onto every tile (except DAPI) when it is saved and 
        % removed from every tile when loaded so we can have negative pixel 
        % values. Saves as uint16.
        TilePixelValueShift = 15000;        
        
        %Tophat filtering is done to each tile in the Dapi channel
        %with a filter of radius DapiR. If set to 'auto', DapiR =
        %round(8/pixelsize). R should be about radius of object of
        %interest.
        DapiR = 'auto';
        
        %MaxWaitTime is the maximum time in seconds that is waited to see
        %if the input data is obtained on the fly. MaxWaitTime1 is for
        %round 1, default value is less so don't wait ages if name is wrong
        %etc.
        MaxWaitTime1 = 60;
        MaxWaitTime = 21600;
        
        %if nPixelsOutsideTiffRange(t,c,r)>nPixelsOutsideTiffRangeThresh
        %then an error will thrown in extract_and_filter.
        nPixelsOutsideTiffRangeThresh = 5000;
        
        %nPixelsOutsideTiffRange(t,c,r) is the number of pixels in tile t,
        %colour channel c, round r that are above the Tiff limit set by
        %uint16(inf). Ideally want all of these to be 0.
        nPixelsOutsideTiffRange;
        
        %PixelsOutsideTiffRangeExtractScale(t,c,r) is the max value that
        %o.ExtractScale can be for all pixels in tile t, round r, colour
        %channel c to be within the Tiff range. If
        %nPixelsOutsideTiffRange(t,c,r) = 0, this is nan
        PixelsOutsideTiffRangeExtractScale;

        %When converting the 3D data to 2D using fstack_modified, all z
        %planes above this will be used.
        FirstZPlane = 2;        
        
        %if StripHack is true, will correct for vertical strips in raw data
        %that are all the same value. In the filtered image, these strips
        %will have a value of 0. 
        StripHack = true;
        
        %% parameters: registration and alignment
        
        % correlation threshold for image alignment. Can be low since 
        % chance correls are miniscule. Second value is for shifts of 0 
        % which can occur spuriously due to hot pixels so need a stronger
        % threshold
        RegCorrThresh = [.3 .6]; 
        
        % smoothing before registration (size of disk filter)
        RegSmooth = 1;
                
        % minimum size overlap for tile matching (pixels - so side squared)
        RegMinSize = 200^2; 
        
        % distance scale for point cloud registration (pixels)
        PcDist = 3; 
        
        %PcIter is the max number of iterations done by the PCR algorithm.
        PcIter = 100;
        
        %if fraction of images in PCR, whose neighbours have converged is
        %less than PcCovergedImgFrac, a warning is given.
        PcCovergedImgFrac = 0.1;
        
        % fraction of point cloud matches needed to count an overlap
        MinPCMatchFract = 0.1; 
        MaxPCMatchThresh = 500;    %required matches can't be greater than this.
        
        % If the number of images (Total number of images =
        % nTiles*nBP*nRounds) with nMatches < MinPCMatches exceeds
        % PcImageMatchesThresh then an error is thrown
        % If auto, set to 1.5*nTiles
        PcImageMatchesThresh = 'auto';
        
        % ToPlot: [r,c,t], plot of round r, colour channel c, tile t
        % is shown for debugging
        ToPlot;
        
        % InitialShiftChannel is the chosen channel to use to register all
        % rounds to the anchor. Choose the one with clear spots, recommend 5,6 or 7.
        InitialShiftChannel = 4;
        
        % MaxRoundShift is the maximum absolute shift in either direction
        % of a tile in any round relative to the anchor round.
        MaxRoundShift = 100;       
        
        % Registration is forced to find an overlap smaller than MaxOverlapFract*TileSz
        % between neighbouring Tiles
        MaxOverlapFract = 0.2;
        
        % MaxRegShift is the maximum shift of a tile in the non overlapping
        % direction i.e. South/North if looking at its east neighbour.
        MaxRegShift = 50;
        
        %Below are required if using register2.m and find_spots2.m - they
        %are to do with the initial search.
        
        %The Score used to find the best shift in get_initial_shift2 is 
        %sum(exp(-Dist.^2/(2*o.ShiftScoreThresh^2)))
        ShiftScoreThresh = 2;
        
        % If either RegMinScore or FindSpotsMinScore are set to 'auto', 
        % the minimum allowed score before the search range is enlarged
        % will be set to median + InitalShiftAutoMinScoreParam*IQR
        % of search scores.
        InitalShiftAutoMinScoreParam = 5;
        
        %RegSearch.Direction.Y,RegSearch.Direction.X are the ranges values of
        %shifts to check during the registration to the neighbour in the corresponding
        %Direction (South or East)
        RegSearch;
        
        %RegStep specifies the step size to use in the Y,X direction of RegSearch
        RegStep = [5,5];               
        
        %if the score mentioned above is below RegMinScore, the search
        %range will be enlarged. If set to 'auto', will be set to median + 5*IQR
        %of search scores.
        RegMinScore = 'auto';
        
        %RegWidenSearch specifies how much to widen the search range in the
        %Y,X directions respectively if the score is below MinRegScore.
        %The search will go from min(RegSearch):RegStep:max(RegSearch) to
        %min(RegSearch)-RegWidenSearch:RegStep:max(RegSearch)+RegWidenSearch
        RegWidenSearch = [50,50];
        
        %After the initial broad search, get_initial_shift will do a
        %refined search with a range of
        %BestShift-RegRefinedSearch:BestShift+RegRefinedSearch inclusive.
        RegRefinedSearch = [12,12];
        
        %RegisterRefinedStep is the step size to use in this refined
        %search
        RegisterRefinedStep = [3,3];
        
        %if the score is below RegAbsoluteMinScore, the shift found will be
        %set to the average of all the other shifts
        RegAbsoluteMinScore = 4;
        
        %To be considered an outlier, a shift must have a score less than
        %OutlierMinScore. AmendShifts will not run unless atleast one of
        %the shifts has a score less than this.
        OutlierMinScore = 100;
        
        %OutlierThresh is the number of scaled MAD away from the median
        %that constitutes an outlier when considering the shifts in the
        %registration or find_steps methods. These outliers are set to the
        %average of all acceptable shifts if the score is low enough
        OutlierThresh = 5;
        
        %RegMaxBadShifts is the number of shifts outside range defined by
        %RegSearch that is considered unacceptable and will cause the
        %method to change to Fft.
        RegMaxBadShifts = 10;
        
        %if the fraction of tiles left to complete is less than 
        %RegInSignificantFract, then the registration step won't start
        %again using the Fft method.
        RegInsignificantFract = 0.1;
        
        %RegMethod is 'PointBased' or 'Fft', if 'PointBased', will only use
        %points detected, if 'Fft' will use convolution method on full
        %images.
        RegMethod = 'PointBased';                
        
        %% parameters: spot detection
        
        % smooth images before reading out fluorescence and before registration 
        % with a disk of this radius. Set to 0 to do no smoothing
        SmoothSize = 0;
        
        % to detect spot, pixel needs to be above dilation with this radius
        DetectionRadius = 1;
        
        % and pixel raw fluorescence needs to be above this value:
        DetectionThresh = 'auto';
        
        % AutoThresh(t,c,r) is the threshold found automatically for tile
        % t, color channel c, round r. Its value is
        % AutoThreshMultiplier*Median(abs(ScaledFilteredTile)).
        AutoThresh;
        AutoThreshMultiplier = 10;
        
        %If find less than o.minPeaks spots then DectionThresh is lowered by
        %o.ThreshParam        
        minPeaks = 1000;
        ThreshParam = 5;
        
        % find isolated spots by annular filtering with these radii
        IsolationRadius1 = 2;
        IsolationRadius2 = 7;
        
        % annular filtered value needs to be less than this:
        IsolationThresh = 'auto';
        
        % for visualization during spot detection
        FindSpotsRoi = [1742 1755 213 227];
        
        % Each spot will be allocated to home tile if possible - but not if
        % it is this close to the edge, because of aberrations
        ExpectedAberration = 3;
        
        %MinThresh is the smallest value DectionThresh can get to       
        MinThresh = 10;
        
        %MinSpots is the smallest number of spots on a round/tile for a
        %particular colour channel for hat color channel to be deemed
        %suitable for finding the initial shifts to the anchor round. If
        %all color channels have a tile with less spots than this, an error
        %is thrown.
        %Also any transforms that have nMatches<MinSpots will be
        %recalculated in AmendPCR
        MinSpots = 50;
        
        %FindSpotsSearch.Y,FindSpotsSearch.X are the ranges values of
        %shifts to check when looking for the initial shifts between rounds
        %for each tile. Can also set to cell(o.nRounds,1) and give a
        %different search range for each round.
        FindSpotsSearch;
        
        %FindSpotsStep specifies the step size to use in the Y,X direction of FindSpotsSearch
        FindSpotsStep = [5,5];
        
        %if the score mentioned above is below FindSpotsMinScore, the search
        %range will be enlarged. If set to 'auto', will be set to median + 5*IQR
        %of search scores.
        FindSpotsMinScore = 'auto';
        
        %RegWidenSearch specifies how much to widen the search range in the
        %Y,X directions respectively if the score is below MinRegScore.
        %The search will go from min(FindSpotsSearch):RegStep:max(FindSpotsSearch) to
        %min(FindSpotsSearch)-FindSpotsWidenSearch:FindSpotsStep:max(FindSpotsSearch)+FindSpotsWidenSearch
        FindSpotsWidenSearch = [50,50];
        
        %After the initial broad search, get_initial_shift will do a
        %refined search with a range of
        %BestShift-FindSpotsRefinedSearch:BestShift+FindSpotsRefinedSearch
        %inclusive.
        FindSpotsRefinedSearch = [12,12];
        
        %FindSpotsRefinedStep is the step size to use in this refined
        %search
        FindSpotsRefinedStep = [3,3];
        
        %if the score is below FindSpotsAbsoluteMinScore, the shift found will be
        %set to the average of all the other shifts in a particular round
        FindSpotsAbsoluteMinScore = 4;
        
        %FindSpotsMethod is 'PointBased' or 'Fft', if 'PointBased', will only use
        %points detected, if 'Fft' will use convolution method on full
        %images to find initial shifts between tiles.
        FindSpotsMethod = 'PointBased';
        
        % PointCloudMethod=1 finds a transform to each tile,round,channel.
        % PointCloudMethod=2 finds a transform to each tile and round and
        % then a scaling to each colour channel. 
        % Use 2 if no anchor round else use 1.
        PointCloudMethod = 1;
        
        
        %% parameters: spot calling
        % normalizes spot fluorescence so this percentile = 1
        SpotNormPrctile = 98;
        
        % if BleedMatrixType == 'Separate', then a bleed matrix will be
        % computed for each round. If BleedMatrixType == 'Single', a single
        % bleed matrix will be computed, combining spot colours from all
        % rounds.
        BleedMatrixType = 'Single';
        
        % This controls how to normalise the spot and bled codes in o.call_spots.
        % If CallSpotsCodeNorm == 'Round', each round of code will have
        % L2 norm of 1. Otherwise, whole code will have L2 norm of 1.
        % Might want to use 'Round', as this means the contribution to each round
        % is the same which is what we expect from the UnbledCodes.
        CallSpotsCodeNorm = 'WholeCode';
        
        % score and intensity thresholds to plot a spot (combi codes)
        %Note max score is 1 for  CallSpotsCodeNorm = 'WholeCode' but 7 for CallSpotsCodeNorm = 'Round'   
        CombiQualThresh = .8;     
        CombiIntensityThresh = .1;
        CombiDevThresh = 0.07;
        CombiAnchorsReq = 4; % need at least this many anchor chans above threshold
        
        nRedundantRounds = 0;
        RedundantPseudobases = {'AC', 'GT'}; 
%         RedundantCodes = {...
%                 '.11..', '.12..', '.22..', '.21..';...
%                 '..11.', '..12.', '..22.', '..21.';...
%                 '1..1.', '2..1.', '2..2.', '1..2.';...
%                 '11...', '12...', '22...', '21...';...
%             };

        RedundantCodes = {...
                '.11..', '.12..', '.21..', '.22..';...
                '..11.', '..12.', '..21.', '..22.';...
                '1..1.', '2..1.', '1..2.', '2..2.';...
                '11...', '12...', '21...', '22...';...
            };

        
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
        
        % for pie-plots: don't both showing probs below this
        MinPieProb = .1;
        
        % size of pies
        PieSize = 25;
        
        % for plotting: any cell classes starting with these names will be
        % collapsed and given the specified color
        ClassCollapse = {{'Astro', 'Endo', 'Oligo', 'Eryth', 'Vsmc'}, 'NonNeuron', [.5 .5 .5] ; ...
                         {'PC.CA1'}, 'PC CA1', [1 .8 .8] ; ...
                         {'PC.CA2', 'PC.CA3'}, 'PC subtypes', [.8 1 .8] ; ...
%                          {'PC.CA3'}, 'PC CA3', [.8 1 .8] ; ...
                         {'Zero'}, 'Zero', [0 0 0]};


        
        %% parameters: stuff that should be the same between experiments
        
        % for dapi images: scale is downsampling factor for final .fig file
        DapiChannel = 1;  
        
        % which channel of each file is anchor images
        AnchorChannel = 2;
        
        % which channel has first round of sequencing images
        FirstBaseChannel = 3;
        
        % which sequencing round to align all others to
        ReferenceRound = 8;
        
        % AnchorRound is index of anchor imaging round that contains Dapi
        % image in DapiChannel.
        AnchorRound;
        
        % ReferenceChannel is the channel in the ReferenceRound that the
        % global coordinate system is built upon.
        ReferenceChannel;
        
        % ReferenceSpotChannels are the channels that are used to find spots
        % for the reference image. If RefRound = Anchor then this is
        % o.AnchorChannel, otherwise it is all channels.
        ReferenceSpotChannels;
        
        % how many combinatorial sequencing rounds excluding anchor round
        nRounds = 7;
        
        % Number of possible basepairs (always 4 for life as we know it but
        % hardcoding is bad programming style!)
        nBP = 7;
        
        % tile size (assumed square tiles)
        TileSz = 2048;
        
        % expected step size (in coordinates returned by bioformats)
        MicroscopeStepSize = 2048;
        
        RawFileExtension = '.nd2';      %e.g. '.nd2', '.czi'
        
        % when decoding spots, will use all colour channels in
        % UseChannels (Array of numbers in range 1 to o.nBP)
        UseChannels;
        
        % when decoding spots, will use all rounds in
        % UseRounds (Array of numbers in range 1 to o.nRounds)
        UseRounds;

        
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
        
        % BigAnchorFile is full path of stitched Anchor Channel image (tiff format)
        BigAnchorFile;
        
        % CellMapFile is full path of cell map (.mat format)
        CellMapFile;
        
        %% variables: registration
        
        % TilePosYX(t,:): grid position of tile t in integers. t is what you
        % find in the file name - so it only counts non-empty tiles. This
        % is confusing and ought to be refactored. Manana.
        TilePosYX;
        
        % TileOriginYX(t,:,r): YX coordinate of origin of tile t on round
        % r. To get the global coordinate of a point, add this to the local
        % coordinate within the tile (counting from 1)
        TileOrigin;
        
        % TileInitialPosYX(t,:): coordinate of tile t in integers. This is
        % the derived coordinates, before using the fact we know the answer
        % to get TilePosYX.
        TileInitialPosYX;
        
        %RawLocalYX{t,b} stores the YX coordinates of spots found in the
        %anchor round of tile t, colour channel b
        RawLocalYX;
        
        %RawIsolated{t,b} labels each spot in the anchor round as isolated or not
        RawIsolated;
        
        %RawLocalNo(t) is the number of spots in tile t of reference round,
        %accounting for all colour channels
        RawLocalNo;
        
        %RegInfo saves debugging information for the registration section
        %h = horizontal, v = vertical
        %Score(i) is the score for the Shifts(i) found between tiles given by Pairs(i)
        %ChangedSearch is number of times the search range had to be changed.
        %Outlier is a shift found that was wrong and subsequently changed
        %to the average of all other shifts.
        RegInfo;      
        
        %AllBaseSpotNo(t,c,r) is the number of spots found on tile t,
        %color channel c, round r.
        AllBaseSpotNo;
        
        % D0(t,2,r) stores the initial shift to use as a starting point for
        % the PCR on round r tile t.
        D0;
        
        %FindSpotsInfo saves debugging information for finding the initial
        %shifts. 
        %Score(t,r) is the score for the shift D0(t,r) found for tile t 
        %between the anchor round and round r.
        %ChangedSearch is number of times the search range had to be changed.
        %Outlier is a shift found that was wrong and subsequently changed
        %to the average of all other shifts.
        %DFailed is D before AmendPCR.
        %nMatchesFailed is nMatches before AmendPCR.
        %ErrorFailed is Error before AmendPCR.
        FindSpotsInfo;
        
        %TileCentre is the yx centre of tile, required for transformation.
        TileCentre;
        
        % A(c): stores the mean scaling correction for chromatic aberration
        % found by point cloud registration for color channel c. This is
        % only used for those for which PcFailed(t,c,r)==true.
        A;
        
        % D(3,2,t,r,c): stores the final affine transform found by point cloud registration
        % on tile t, round r, channel c.
        D;
        
        % cc(t,1,r) stores the correlation coefficient for the initial
        % shift D0(t,:,r) found between tile t round r and the anchor
        cc;
        
        % nMatches(t,c,r): stores number of matches found by point cloud
        % registration for tile t, color channel c, round r
        nMatches;
        
        % error(t,c,r): stores error found by point cloud registration
        % for tile t, color channel c, round r
        Error;
        
        % PcFailed(t,c,r): This is a logical array. A value of 1 indicates
        % that nMatches high enough for the PCR results to be used. For
        % these images, chromatic aberration will be the mean value of that
        % colour channel. The shift will be found using an Fft method.
        PcFailed;
                
        % nPcCovergedImg is the fraction of images that converged in PCR.
        % Denominator is nTiles*o.nBP*o.nRounds
        nPcCovergedImg;
        
        % PcRegLambda is the constant for the regularised PCR. 1 is the
        % shift scale factor. 2 is the scaling (between colour channels)
        % scale factor
        PcRegLambda1 = 9;
        PcRegLambda2 = 30000;
        
        % PcMinSpots is the minimum number of spots to do PCR without
        % regularization. Below this amount regularisation is used.
        PcMinSpots = 200;       
        
        % If a tile has a chromatic aberration scaling that has an 
        % absolute deviation of more than PcMaxScaleDev from the median 
        % for that colour channel and nMatches<PcMinSpotsScaling,
        % it will be re-evaluated with regularization.
        PcMaxScaleDev = 0.5e-3;
        PcMinSpotsScaling = 500;
        
        % If a tile has a shift that has an absolute deviation of more than
        % PcMaxShiftDev pixels from the median for that tile and round, it will be
        % re-evaluated with regularization.
        PcMaxShiftDev = 10;
        
        % PcGrad(t,c,r,1) is the gradient of the Y shift with respect to
        % the normalised Y coordinate (from -1 to +1) for tile t, channel
        % c, round r. PcGrad(t,c,r,2) is the X gradient with respect to x
        % coordinate. abs value should be <0.5 if worked.
        PcGrad;
        
        % PcMean(t,c,r,1) is the mean of the Y shift with respect to
        % the normalised Y coordinate (from -1 to +1) for tile t, channel
        % c, round r. PcGrad(t,c,r,2) is the X mean with respect to x
        % coordinate. abs value should be <0.5 if worked.
        PcMean;

        %% variables: dot product spot calling outputs
        
        % dpSpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        dpSpotGlobalYX;
       
        % dpSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        dpSpotColors;             
        
        % dpSpotIsolated(Spot) is a binary array saying if the spot is well isolated
        % again for combinatorial spots only
        dpSpotIsolated; 
        
        % dpLocalTile(s) is the tile spot cSpotColors(s,:,:) was found on
        dpLocalTile;
        
        % SpotCombi: binary array which is one for combinatorial spots,
        % zero for extra spots
        dpSpotCombi;       
        
        % dpSpotCodeNo(Spot): code number for each spot (uint16). Both combinatorial and extra spots
        dpSpotCodeNo;        
     
        % dpSpotScore(Spot): score saying how well the code fits (0...1).
        % 1 for extras
        dpSpotScore;
        
        % dpSpotIntensity(Spot): RMS intensity of the spot. Zero for
        % extras!
        dpSpotIntensity;
        
        % dpSpotScoreDev(s) is the standard deviation of the scores of
        % assigning spot s to all genes
        dpSpotScoreDev;
        
        % CharCodes{Code}: text representation of each code in the
        % codebook. For extra spots, this just says "EXTRA"
        CharCodes;
        
        % GeneNames{Code}: gene name for each code. Note that some genes
        % can have multiple entries, because they have multiple codes or
        % because they have extra rounds
        GeneNames;
        
        % UnbledCodes(nCodes, nBP*nRounds): binary code vectors
        UnbledCodes;
        
        % BledCodes(nCodes, nBP*nRounds): code vectors after modeling
        % crosstalk
        BledCodes;
        
        % Normalised Spot Scores
        NormBledCodes;
        dpNormSpotColors;
        
        % BleedMatrix used to estimate BledCodes
        BleedMatrix;
        
        % SpotColors./BledCodesPercentile are the spot colors used to get
        % BleedMatrix. 
        BledCodesPercentile;
		
        %% variables: cell calling outputs
        % pCellClass(cell, class); % prob each cell goes to each class: last class is zero expression
        pCellClass;
        
        % ClassNames(class): name of each class, taken from gSet file
        ClassNames;
        
        % pSpotCell(spot, cell): sparse array containing prob of each spot
        % to belong to each cell. Last cell is background
        pSpotCell;
        
        % position of each cell centroid
        CellYX;
        
        
        %% parameters: Bckground distribution for Probability Method of spot calling        
        %HistCounts(:,b,r) is the pixel count corresponding to each value
        %in HistValues for channel b, round r
        HistCounts;
        
        %HistValues: The full range of pixel values;
        HistValues;
        
        %HistMaxValues(b,r) is the pixel value where the histogram for
        %channel b, round r is peaked. Expect it to be 0.
        HistMaxValues;
        
        %if the maximum absolute value of HistMaxValues is greater than
        %HistMaxShiftThresh then an error is thrown.
        HistMaxShiftThresh = 25; 
        
        %% Prefixes for the different spot calling methods
        CallMethods = {'DotProduct','Prob','Pixel','OMP','Spatial','GroundTruth'};
        % DotProduct: Takes all spots on anchor round, score is dot product
        % between normalized spot colors and bled codes.
        % Prob: Takes all spots on anchor round, score for each
        % round/channel is probability it can be explained by gene and
        % background minus probability background only can explain it.
        % Pixel: Same as prob but look at all pixels not just anchor spots.
        % OMP: for each pixel, iteratively fits genes until extra genes
        % don't reduce the residual by more than certain amount.
        % Spatial: Like OMP, but now fitting genes with a spatial
        % distribution not just singular pixels. 
        CallMethodPrefix = ...
            containers.Map({'DotProduct','Prob','Pixel','OMP','Spatial','GroundTruth'},...
            {'dp','p','px','omp','sp','gt_px'});
    end
    
end

 





