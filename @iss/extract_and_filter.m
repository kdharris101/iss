function o = extract_and_filter(o)
% create tiff files for each tile that are top-hat filtered versions of
% original czi files
    GPU_test = gpuArray([1]);       %So if no Parallel Computing Toolbox, fails straight away
    o.TileFiles = cell(o.nRounds+o.nExtraRounds,1,1); % 1,1 because we don't yet know how many tiles
    
    for r = 1:o.nRounds+o.nExtraRounds       
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);

        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        
        if exist(imfile)>0
            % initiate reader
            bfreader.setId(imfile);
        else
            %Wait for file to exist
            if r==1
                MaxTime = o.MaxWaitTime1;   %Don't wait long if first round
            else
                MaxTime = o.MaxWaitTime;   %Max amount of time to wait in seconds
            end
            count = 0;
            while exist(imfile)==0
                pause(1);
                count = count+1;
                if count >= MaxTime
                    error(sprintf(strcat('No file named:\n  ',imfile,'\ncreated in allowed time')));
                end
            end
            %Wait for file to stop loading
            OldBytes = 0;
            NewBytes = 0.00001;
            while NewBytes>OldBytes
                pause(5);
                fileinfo = dir(imfile);
                OldBytes = NewBytes;
                NewBytes = fileinfo.bytes;
            end
            % initiate reader
            bfreader.setId(imfile);
        end

        % get some basic image metadata
        [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
            get_ome_tilepos(bfreader);
        if isempty(xypos) || size(xypos, 1)==1
            if r == 1
                warning('first round xypos empty - using values from initial manual input')
                assert(~isempty(o.TileInitialPosYX), 'xypos unavailable')
                xypos = o.TileInitialPosXY;
                xyposOld = xypos;
            else
                warning('xypos empty - using values from previous round')
                xypos = xyposOld;
            end
            nSerieswPos = size(xypos,1);
        else
            xyposOld = xypos;
        end
        
        scene = nSeries/nSerieswPos;

        bfreader.close();
        
        if r == 1
            if isempty(o.AutoThresh)
                o.AutoThresh = zeros(nSerieswPos,nChannels,o.nRounds+o.nExtraRounds);
            end
            if isempty(o.HistCounts)
                o.HistValues = -o.TilePixelValueShift:1:2^16-o.TilePixelValueShift;  %Entire range of possible pixel values
                o.HistCounts = zeros(length(o.HistValues),nChannels,o.nRounds);
            end
            if isempty(o.TilePosYX)
                % find x and y grid spacing as median of distances that are about
                % right
                dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
                xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
                dy = xypos(:,1)-xypos(:,1)'; % all pairs of y distances
                yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));
        
        
                % find coordinates for each tile
                if isempty(o.TileInitialPosYX)
                    o.TileInitialPosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));
                end
                o.TilePosYX = o.TileInitialPosYX;
                %Below is a safeguard incase wrong positions found - can do
                %this as we knwo what the answer should be.
                MaxY = max(o.TileInitialPosYX(:,1));
                MaxX = max(o.TileInitialPosYX(:,2));
                if MaxY*MaxX ~= nSeries
                    warning('Number of tiles (%d) is not equal to maximum Y position (%d) multiplied by maximum X position (%d)'...
                        , nSeries, MaxY, MaxX)
                    break
                else
                    TilePosY = flip(repelem(1:MaxY,MaxX));
                    o.TilePosYX(:,1) = TilePosY;
                    TilePosX = repmat([flip(1:MaxX),1:MaxX],1,ceil(MaxY/2));
                    o.TilePosYX(1:nSeries,2) = TilePosX(1:nSeries);
                end
            end
            
            %New filter
            if strcmpi(o.ExtractR1, 'auto') || strcmpi(o.ExtractR2, 'auto')
                o.ExtractR1 = round(0.5/pixelsize);  %Gives value of 3 for pixelsize = 0.1669 of most data tested
                o.ExtractR2 = o.ExtractR1*2;
            end
            h = -hanning(o.ExtractR2*2+1);
            h = -h/sum(h);
            h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1) = ...
                h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1)+hanning(o.ExtractR1*2+1)/sum(hanning(o.ExtractR1*2+1));
            SE = ftrans2(h');
            SE = single(gpuArray(SE));
            
        end

        % set up filename grid for this round
        fName = cell(nSerieswPos,1);
        
        %Set top hat structuring elements
        if strcmpi(o.DapiR,'auto')
            o.DapiR = round(8/pixelsize);
        end
        DapiSE = strel('disk', o.DapiR);
        
        %parfor t = 1:nSerieswPos  
        for t = 1:nSerieswPos  
           
            fName{t} = fullfile(o.TileDirectory, ...
                    [o.FileBase{r}, '_t', num2str(t), '.tif']);
            
            if exist(fName{t}, 'file')
                fprintf('Round %d tile %d already done.\n', r, t);
                if o.AutoThresh(t,o.AnchorChannel,r) == 0
                    TifObj = Tiff(fName{t});
                    for c=1:nChannels
                        if c ~= o.AnchorChannel && r == o.ReferenceRound; continue; end
                        TifObj.setDirectory(o.FirstBaseChannel + c - 1);
                        IFS = int16(TifObj.read())-o.TilePixelValueShift;
                        o.AutoThresh(t,c,r) = median(abs(IFS(:)))*o.AutoThreshMultiplier;
                        if ismember(r,1:o.nRounds)
                            o.HistCounts(:,c,r) = o.HistCounts(:,c,r)+histc(IFS(:),o.HistValues);
                        end
                    end
                end
                continue;
            end                   
            
                
            % a new reader per worker
            bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % use the memo file cached before
            bfreader.setId(imfile);

            bfreader.setSeries(scene*t-1);

            for c = 1:nChannels
                % read z stacks
                tic
                I = cell(nZstacks,1);
                for z = 1:nZstacks
                    iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                    I{z} = bfGetPlane(bfreader, iPlane);
                end

                % focus stacking
                I_mod = o.fstack_modified(I);

                % tophat
                if c == o.DapiChannel && r == o.ReferenceRound
                    IFS = imtophat(I_mod, DapiSE);
                else
                    I_mod = single(padarray(I_mod,(size(SE)-1)/2,'replicate','both'));
                    IFS = convn(gpuArray(I_mod),SE,'valid'); 
                    clearvars I_mod I  %Free up GPU memory
                    
                    if strcmpi(o.ExtractScale, 'auto')
                        o.ExtractScale = round(30000/max(IFS(:)));
                    end
                    IFS = IFS*o.ExtractScale;
                    
                    if c ~= o.AnchorChannel && r == o.ReferenceRound
                        IFS = gather(uint16(IFS+o.TilePixelValueShift)); 
                    else
                        %Determine auto thresholds
                        o.AutoThresh(t,c,r) = gather(median(abs(IFS(:)))*o.AutoThreshMultiplier);
                        if ismember(r,1:o.nRounds)
                            %Get histogram data
                            IFS = int16(IFS);
                            o.HistCounts(:,c,r) = o.HistCounts(:,c,r)+gather(histc(IFS(:),o.HistValues));
                        end
                        IFS = gather(uint16(IFS+o.TilePixelValueShift)); 
                    end

                end

                % write stack image
                imwrite(IFS,...
                    fullfile(o.TileDirectory,...
                    [o.FileBase{r}, '_t', num2str(t), '.tif']),...
                    'tiff', 'writemode', 'append');
                toc
            end
            fprintf('Round %d tile %d finished.\n', r, t);
            bfreader.close();

        end        

        for t=1:nSerieswPos
            o.TileFiles{r,o.TilePosYX(t,1), o.TilePosYX(t,2)} = fName{t};
        end
    end
    
    o.EmptyTiles = cellfun(@isempty, squeeze(o.TileFiles(o.ReferenceRound,:,:)));

    %Plot boxplots showing distribution af AutoThresholds
    if o.Graphics
        UseRounds = setdiff(1:o.nRounds+o.nExtraRounds,o.ReferenceRound);
        Thresholds = [];
        group = [];
        index = 1;
        for c=1:nChannels
            for r=UseRounds
                Thresholds = [Thresholds;o.AutoThresh(:,c,r)];
                group = [group;index*ones(size(o.AutoThresh(:,1,1)))];
                index = index+1;
            end
        end
        %Add anchor
        AnchorLabel = {'Anchor'};
        Thresholds = [Thresholds;o.AutoThresh(:,o.AnchorChannel,o.ReferenceRound)];
        group = [group;index*ones(size(o.AutoThresh(:,1,1)))];

        
        figure(43290);
        colors = colormap(lines(nChannels));
        Colors = repelem(colors,length(UseRounds),1);
        Colors = [Colors;repelem([0,0,0],nChannels,1)];
        Labels = [string(repmat(UseRounds,1,nChannels)),string(AnchorLabel)];
        boxplot(Thresholds,group,'Colors',Colors, 'plotstyle', 'compact','labels', Labels);
        set(gca,'TickLength',[0 0]);
        ylabel('AutoThreshold');
        xlabel('Round');
        hold on
        for c=1:nChannels
            plot(NaN,1,'color', colors(c,:), 'LineWidth', 4);       %For legend labels
        end
        leg = legend(o.bpLabels,'Location','northwest');
        title(leg,'Color Channel');
        hold off
    end
    
    %Plot histograms to make sure they are smooth
    %Avoid ExtraRounds as only need histograms for the 7 rounds used to
    %define genes
    if o.Graphics       
        figure(43291);
        index = 1;
        for r=1:o.nRounds
            for b=1:nChannels
                subplot(o.nRounds,nChannels,index)
                histogram('BinEdges',[o.HistValues-0.5,max(o.HistValues)+0.5],'BinCounts',o.HistCounts(:,b,r),'DisplayStyle','stairs');
                xlim([-1000,1000]);
                ylim([0,max(o.HistCounts(:,b,r))]);
                if b==4
                    title(strcat('Round ',num2str(r)));
                end
                index = index+1;
            end
        end
    end
end
