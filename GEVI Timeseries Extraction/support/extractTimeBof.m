function timeBof = extractTimeBof(filename)
    metadata = parseTiffHeader(filename);
    metadataLines = strsplit(metadata, '\n');
    if length(metadataLines)>20 && contains(metadataLines{21}, 'meta.header.timeBof')
        timeBof = str2double(regexp(metadataLines{21}, '\d+', 'match'));
        return
    end
    for i = 1:length(metadataLines)
        if contains(metadataLines{i}, 'meta.header.timeBof')
            timeBof = str2double(regexp(metadataLines{i}, '\d+', 'match'));
            return;
        end
    end
    error('meta.header.timeBof not found in file %s', filename);
end


function imageDescription = parseTiffHeader(filePath)
    % Open the TIFF file
    fid = fopen(filePath, 'r'); % 'b' for big-endian format

    % Check if file is opened successfully
    if fid == -1
        error(['Could not open file: ', filePath]);
    end

    % Read Byte Order (II or MM)
    byteOrder = fread(fid, 2, 'char')';

    if isequal(byteOrder, 'II')
        endianString = 'l';
    elseif isequal(byteOrder, 'MM')
        fclose(fid);
        error('Big-endian TIFF files are not supported.');
    else
        fclose(fid);
        error('Invalid TIFF byte order');
    end

    % Read the next 2 bytes (magic number) as a uint16
    fseek(fid, 2, 'bof'); % Move to the beginning of the file
    magicNumber = fread(fid, 1, 'uint16',0,endianString);

    if magicNumber ~= 42
        fclose(fid);
        error('Invalid TIFF magic number');
    end

    % Read the first IFD offset (4 bytes)
    ifdOffset = fread(fid, 1, 'uint32', 'l');

    % Close the file
    fclose(fid);

    % Read Image Description
    fid = fopen(filePath, 'r', 'l'); % Reopen in little-endian format

    if fid == -1
        error(['Could not open file: ', filePath]);
    end

    % Move to the IFD offset
    fseek(fid, ifdOffset, 'bof');

    % Read the number of entries in the IFD
    numEntries = fread(fid, 1, 'uint16', 0, 'l');

    % Iterate through IFD entries
     for i = 1:numEntries
        tag = fread(fid, 1, 'uint16', 0, 'l');
        if tag == 270 % Image Description tag
            fieldType = fread(fid, 1, 'uint16', 0, 'l');
            numValues = fread(fid, 1, 'uint32', 0, 'l');
            valueOffset = fread(fid, 1, 'uint32', 0, 'l');

            % Move to the position of the value
            fseek(fid, valueOffset, 'bof');

            % Read the Image Description
            imageDescription = fread(fid, numValues, '*char')';

            break; % No need to continue searching
        else
            % Skip to the next IFD entry
            fseek(fid, 10, 'cof'); % Each IFD entry is 12 bytes
        end
    end

    % Close the file
    fclose(fid);
end
