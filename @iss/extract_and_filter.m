function o = extract_and_filter(o)
% create tiff files for each tile that are top-hat filtered versions of
% original czi files

    o.TileFiles = cell(o.nRounds,1,1); % 1,1 because we don't yet know how many tiles

    for r = 1:o.nRounds+o.nExtraRounds
        imfile = fullfile(o.InputDirectory, [o.FileBase{r}, '.czi']);

        % construct a Bio-Formats reader with the Memoizer wrapper
        bfreader = loci.formats.Memoizer(bfGetReader(), 0);
        % initiate reader
        bfreader.setId(imfile);

        % get some basic image metadata
        [nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
            get_ome_tilepos(bfreader);
        if isempty(xypos) || size(xypos, 1)==1
            if r == 1
                warning('first round xypos empty - using values from initial manual input')
                assert(~isempty(o.TileInitialPosXY), 'xypos unavailable')
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
        
        % find x and y grid spacing as median of distances that are about
        % right
        dx = xypos(:,1)-xypos(:,1)'; % all pairs of x distances
        xStep = median(dx(abs(1- dx(:)/o.MicroscopeStepSize)<.5));
        dy = xypos(:,1)-xypos(:,1)'; % all pairs of y distances
        yStep = median(dy(abs(1- dy(:)/o.MicroscopeStepSize)<.5));
        
        % find coordinates for each tile
        o.TilePosYX = fliplr(1+round((xypos - min(xypos))./[xStep yStep]));

        % set up filename grid for this round
        fName = cell(nSerieswPos,1);
        
        parfor t = 1:nSerieswPos  
           
            fName{t} = fullfile(o.TileDirectory, ...
                    [o.FileBase{r}, '_t', num2str(t), '.tif']);
            
            if exist(fName{t}, 'file')
                fprintf('Round %d tile %d already done.\n', r, t);
                continue;
            end                   
            
                
            % a new reader per worker
            bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % use the memo file cached before
            bfreader.setId(imfile);

            bfreader.setSeries(scene*t-1);

            for c = 1:nChannels
                % structuring element for top-hat
                if c == o.DapiChannel
                    SE = strel('disk', round(8/pixelsize));     % DAPI
                else
                    SE = strel('disk', round(1/pixelsize));
                end


                % read z stacks
                I = cell(nZstacks,1);
                for z = 1:nZstacks
                    iPlane = bfreader.getIndex(z-1, c-1, 0)+1;
                    I{z} = bfGetPlane(bfreader, iPlane);
                end

                % focus stacking
                IFS = o.fstack_modified(I);

                % tophat
                IFS = imtophat(IFS, SE);

                % write stack image
                imwrite(IFS,...
                    fullfile(o.TileDirectory,...
                    [o.FileBase{r}, '_t', num2str(t), '.tif']),...
                    'tiff', 'writemode', 'append');
            end
            fprintf('Round %d tile %d finished.\n', r, t);
            bfreader.close();

        end        

        for t=1:nSerieswPos
            o.TileFiles{r,o.TilePosYX(t,1), o.TilePosYX(t,2)} = fName{t};
        end
    end
    
    o.EmptyTiles = cellfun(@isempty, squeeze(o.TileFiles(o.ReferenceRound,:,:)));

end
