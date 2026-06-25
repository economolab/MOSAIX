classdef FrameReader < handle
    properties
        fileID
        w
        h
        g
        firstOffset = 0
        dtype
        maxFrames
        frame_size
        shuffle
        pointer
        isTiff
        metadata
        region = [] %currently only for raw files (should be fixed but this is for my use so sucks).
    end
    
    methods
        function obj = FrameReader(fileID, varargin)
            % Initialize input parser
            p = inputParser;
            addRequired(p, 'fileID', @(x) ischar(x) || isstring(x));
        
            % Define default parameters 
            defaultParameters = {'w', 588; 'h', 624; 'g', 1728; 'dtype', 'uint8'; 
                                 'maxFrames', 24000; 'shuffle', false};
                             
            for i = 1:size(defaultParameters, 1)
                addParameter(p, defaultParameters{i, 1}, defaultParameters{i, 2});
            end
        
            % Parse inputs
            parse(p, fileID, varargin{:});
        
            % Assign properties from parsed inputs
            obj.fileID = p.Results.fileID; 
            obj.isTiff = endsWith(obj.fileID, '.tif') || endsWith(obj.fileID, '.tiff'); 
            obj.shuffle = p.Results.shuffle; 
            obj.pointer = 1-2*obj.shuffle; 
            obj.metadata = [];
            % If Tiff, use metadata; if not, use passed or default parameters
            if obj.isTiff
                folder = fileparts(fileID);
                tiffInfo = Tiff(obj.fileID, 'r');
                obj.w = tiffInfo.getTag('ImageWidth'); 
                obj.h = tiffInfo.getTag('ImageLength'); 
                obj.dtype = class(tiffInfo.read());
                try
                    obj.maxFrames = getMetaVal(tiffInfo.getTag('ImageDescription'),'images');
                catch
                    while ~tiffInfo.lastDirectory()
                        tiffInfo.nextDirectory();
                    end
                    obj.maxFrames = tiffInfo.currentDirectory();
                end
                if obj.maxFrames ==-1
                    obj.maxFrames = length(imfinfo(obj.fileID)); 
                end
                obj.g = 0; 
                tiffInfo.close();

            elseif isfolder(obj.fileID)
                folder = obj.fileID;
                % Get a list of all TIFF files in the folder
                
                fileNames = strtrim(string(ls(fullfile(obj.fileID,'*.tiff'))));
                frames = regexp(fileNames,'(\d+)\.tiff','tokens');
                frames = [frames{:}];
                frames = arrayfun(@str2double,[frames{:}]);
                frameLims = [min(frames),max(frames)];
                if ~length(fileNames) == frameLims(2)-frameLims(1)+1
                    tifFiles = dir(fullfile(obj.fileID, '*.tif'));
                    tiffFiles = dir(fullfile(obj.fileID, '*.tiff'));
                    allFiles = [tifFiles; tiffFiles];
                    
                    % Define a function to get imfinfo for a file
                    getImfinfo = @(file) struct(...
                        'fileName', string(fullfile(file.folder,file.name)), ...
                        'tm', getMetaVal(imfinfo(fullfile(obj.fileID, file.name)).ImageDescription,'meta.header.timeBof') ...
                    );
                    
                    % Apply the function to all files
                    imfinfoData = arrayfun(getImfinfo, allFiles);
    %                 if any([imfinfoData.tm]<=0)
                    
                    [~,idx] = sort([imfinfoData.tm]);
                    min_tm = min([imfinfoData.tm]);
                    imfinfoData = arrayfun(@(x) setfield(x, 'tm', x.tm - min_tm), imfinfoData);
    
                    obj.metadata = imfinfoData(idx);
                else
                    fileStart = regexp(fileNames(1), '^(.*?)\d+\.tiff$', 'tokens', 'once');
                    fileStart = fileStart{1};
                    fileNamesOrdered = fileStart + string(frameLims(1):frameLims(2)) + ".tiff";
                    filePath = obj.fileID;
                    tmStart = extractTimeBof(fullfile(filePath, fileNamesOrdered(1)));
                    tmEnd = extractTimeBof(fullfile(filePath, fileNamesOrdered(end)));
                    stepSize = (tmEnd-tmStart)/(length(fileNamesOrdered)-1);
                    tm = tmStart:stepSize:tmEnd;
%                     tm = arrayfun(@(x) extractTimeBof(fullfile(filePath, x)), fileNamesOrdered);

                    obj.metadata = struct('fileName',cellstr(fileNamesOrdered),'tm',num2cell(tm));
                end

                clear imfinfoData
                tiffInfo = Tiff(fullfile(obj.fileID,obj.metadata(1).fileName), 'r');
                obj.w = tiffInfo.getTag('ImageWidth'); 
                obj.h = tiffInfo.getTag('ImageLength'); 
                obj.dtype = class(tiffInfo.read());
                obj.maxFrames = length(obj.metadata);
                obj.g = 0; 
                tiffInfo.close();
            else
                % Check for textFiles and matFiles in the folder
                folder = fileparts(fileID);
                textFiles = dir(fullfile(folder, '*import*.txt'));
                
                % Parse info from textFile if available
                if ~isempty(textFiles) && isscalar(textFiles)
                    [w, h, g, firstOffset, dtype] = parseFileInfo(fullfile(folder,textFiles(1).name)); %#ok<ASGLU>
                end
                
                if firstOffset==80
                    obj.region = getROIRaw(obj.fileID);
                end
                
                % Assign the parameters based on their existence
                parameters = ["w", "h", "g", "firstOffset", "maxFrames", "dtype"]; %there should be a better way of doing this
                for parameter = parameters
                    if exist(parameter, "var")
                        obj.(parameter) = eval(parameter);
                    else
                        obj.(parameter) = p.Results.(parameter);
                    end
                end
            end
            csvFiles = dir(fullfile(folder, '*.csv'));
            matFiles = dir(fullfile(folder, '*.mat'));

            if length(csvFiles)>=1 && any(endsWith(csvFiles.name,"metadata.csv"))
                myCsvFile = csvFiles(find(endsWith(csvFiles.name,"metadata.csv"),1));
                frameTimes = table2array(importfile_CSV_frametimes(fullfile(folder,myCsvFile.name)));
                obj.maxFrames = length(frameTimes);
                obj.metadata.tm = frameTimes;
%                 if obj.maxFrames==
%                     obj.metadata.tm = frameTimes;
%                 else
%                     obj.metadata.tm = (0:obj.maxFrames-1)*mean(diff(frameTimes));
%                 end
            elseif ~isempty(matFiles) && isscalar(matFiles) % Load data from matFile if available
                a = load(fullfile(folder,matFiles(1).name));
                try
                    obj.metadata.tm = a.params.tm;
                catch
                    disp("no metadata for " + obj.fileID + " found!")
                end
                clear a;
            end
        
            % Compute frame size
            obj.frame_size = obj.w * obj.h;
        end


        function frames = getFrames(obj, numFrames, givePartial)
            if nargin < 3
                givePartial = true;
            end
            
            if obj.isDone(numFrames, givePartial)
                frames = [];
                return;
            end
            
            numFrames = min(numFrames,obj.maxFrames-obj.pointer+1);

            start_frame = getStartFrame(obj, numFrames);
            if obj.isTiff
                frames = readTiff(obj, numFrames, start_frame);
            elseif isfolder(obj.fileID)
                frames  = readTiffs(obj, numFrames, start_frame);
            else
                frames = readRaw(obj, numFrames, start_frame);
            end
            obj.pointer = obj.pointer + numFrames;
        end

        function done = isDone(obj, numFrames, givePartial)
            if nargin < 3
                givePartial = true;
            end
            
            done = obj.pointer == obj.maxFrames + 1 || (~givePartial && obj.pointer + numFrames > obj.maxFrames);
        end

        function frames = getFramesByIndexes(obj, startDex,endDex)
            if (endDex>obj.maxFrames)
                error('Invalid frame index: %d', endDex);
            end
            
            if obj.isTiff
                frames = readTiff(obj,endDex-startDex+1,startDex);
            elseif isfolder(obj.fileID)
                frames = obj.readTiffs(endDex-startDex+1,startDex);
            else
                frames = readRaw(obj,endDex-startDex+1,startDex);
            end
        end

    end

    methods (Access = private)
        function start_frame = getStartFrame(obj, numFrames)
            if obj.pointer == -1
                start_frame = randi([0, obj.maxFrames - numFrames]);
            else
                start_frame = obj.pointer;
            end
        end

        function frame = readFrameByIndex(obj, fileObj, frameIndex)
            if obj.isTiff
                fileObj.setDirectory(frameIndex + 1);
                frame = fileObj.read();
            else
                pixelBytes = getPixelBytes(obj.dtype);
                offset = obj.firstOffset;
                offset = offset + (obj.frame_size*pixelBytes + obj.g) * frameIndex;
                fseek(fileObj, offset, 'bof');
        
                image_data = fread(fileObj, obj.frame_size, obj.dtype);
                frame = reshape(image_data, obj.w, obj.h)'; 
            end
        end

        function frames = readTiff(obj, numFrames, start_frame)
%             try
%                 frames = tiffreadVolume(obj.fileID,"PixelRegion",{[1 1 inf],[1 1 inf],[start_frame,start_frame+numFrames-1]});
%                 return
%             catch
%             end
            try
                tiffObj = Tiff(obj.fileID, 'r');
                tiffObj.setDirectory(start_frame);
                
                frames = cell(numFrames, 1);
                for i = 1:numFrames-1
                    frames{i} = tiffObj.read();
                    tiffObj.nextDirectory();
                end
                
                frames{numFrames} = tiffObj.read();
                
                tiffObj.close();
                
                frames = cat(3, frames{:});
            catch ME
                warning(sprintf("using tiffreadVolume because the tiff library failed to read the image. Error: %s", ME.message)); %#ok<SPWRN>
                frames = tiffreadVolume(obj.fileID,"PixelRegion",{[1 1 inf],[1 1 inf],[start_frame,start_frame+numFrames-1]});
                % end
            end
        end

        function frames = readTiffs(obj, numFrames, start_frame)
            directory = obj.fileID;
            frames = imread(fullfile(directory,obj.metadata(start_frame).fileName));
            frames = repmat(frames,[1,1,numFrames]);
            relevantMetadata = obj.metadata(start_frame:start_frame+numFrames-1);
            for i = 1:numFrames
                frames(:,:,i) = imread(fullfile(directory,relevantMetadata(i).fileName));
            end
        end
        

         function frames = readRaw(obj, numFrames, start_frame)
            % move file opening to outside of the function and pass fid as parameter
            fid = fopen(obj.fileID, 'rb'); 
            
            % preallocate array for all frames
            frames = zeros(obj.w, obj.h, numFrames, obj.dtype); 

            pixelBytes = getPixelBytes(obj.dtype);
            
            % calculate the number of batches needed
            batchSize = 1000;
            numBatches = ceil(numFrames / batchSize);
            
            for b = 1:numBatches
                % calculate the number of frames in this batch
                if b == numBatches
                    % for the last batch, the number of frames might be less than batchSize
                    numFramesInBatch = numFrames - (b - 1) * batchSize;
                else
                    numFramesInBatch = batchSize;
                end
            
                % calculate the offset for this batch
                offset = obj.firstOffset;
                offset = offset + (obj.frame_size*pixelBytes + obj.g) * ((b - 1) * batchSize + start_frame - 1);
                fseek(fid, offset, 'bof');
            
                % read data for this batch
                total_elements = numFramesInBatch * obj.frame_size;
                total_gap = (numFramesInBatch - 1) * obj.g;
                total_data_size_pixels = total_elements + total_gap/pixelBytes;
                all_data = fread(fid, total_data_size_pixels, ['*' obj.dtype]);
            
                % create a mask to remove gaps
                gap_step_pixels = obj.frame_size + obj.g/pixelBytes; % step size for each gap
                gap_starts = obj.frame_size + 1 : gap_step_pixels : total_data_size_pixels; % starting indices of each gap
                gap_ends = gap_starts + obj.g/pixelBytes - 1; % ending indices of each gap
                mask = true(total_data_size_pixels, 1); 
                for i = 1:length(gap_starts)
                    mask(gap_starts(i):gap_ends(i)) = false;
                end
                if length(mask)>length(all_data)
                    warning("mask longer than data. Likely fewer frames in file than expected"); %something may be messed up (i.e. data acq interupted?)
                    %basically means metadata and number of frames are messed up
                    mask = mask(1:length(all_data));
                end
                all_data = all_data(mask);
            
                % reshape data into frames
                for i = 1:numFramesInBatch
                    startIdx = (i - 1) * obj.frame_size + 1;
                    endIdx = startIdx + obj.frame_size - 1;
                    image_data = all_data(startIdx : endIdx);
                    frameNum = (b - 1) * batchSize + i; % the overall frame number in frames
                    frames(:,:,frameNum) = reshape(image_data, obj.w, obj.h);
                end
            end
            
            frames = permute(frames, [2 1 3]);
            fclose(fid);
        end
    end
end


function [width, height, gap, firstOffset, imageTypeString] = parseFileInfo(fileInfo)
    % Read the contents of the text file
    fid = fopen(fileInfo, 'r');
    fileInfoContents = fread(fid, '*char')';
    fclose(fid);
    
    % Regular expressions for width, height, gap, and image type
    widthRegex = '- Width: (\d+)';
    heightRegex = '- Height: (\d+)';
    gapRegex = '- Gap between images: (\d+)';
    imageTypeRegex = '- Image Type: ''(\d+)-bit ?\w*''';
    offsetRegex = '- Offset to first image: (\d+)';
    
    % Extract the width, height, gap, and image type
    widthTokens = regexp(fileInfoContents, widthRegex, 'tokens', 'once');
    heightTokens = regexp(fileInfoContents, heightRegex, 'tokens', 'once');
    gapTokens = regexp(fileInfoContents, gapRegex, 'tokens', 'once');
    offsetTokens = regexp(fileInfoContents, offsetRegex, 'tokens', 'once');
    imageTypeTokens = regexp(fileInfoContents, imageTypeRegex, 'tokens', 'once');
    
    % Convert the extracted tokens to numeric values
    width = str2double(widthTokens{1});
    height = str2double(heightTokens{1});
    gap = str2double(gapTokens{1});
    firstOffset = str2double(offsetTokens{1});
    imageTypeString = ['uint' imageTypeTokens{1}];
end

function region = getROIRaw(fn)
% Open the file
fid = fopen(fn, 'r', 'l');

% Read the first 80 bytes
data = fread(fid, 80, 'uint16');

% Close the file
fclose(fid);

region = [data(30)+1,data(31)+1,data(33)+1,data(34)+1];
end

function pixelBytes = getPixelBytes(dtype)
    if strcmp(dtype,'uint8')
        pixelBytes = 1;
    elseif strcmp(dtype,'uint16')
        pixelBytes = 2;
    else
        error("doesn't support that datatype")
    end
end
