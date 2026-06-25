function paramsArray = Select_input_files(app, rootDir)
    if ~exist("rootDir",'var')
        rootDir = getpref(class(app),'lastPath',pwd);
    end
    tempParam.framesToDisregard = 60;
    tempParam.Nplanes=1;
    tempParam.isRaw=1;
    tempParam.vertSplit=0;
    
    prompt = 'Select image(s) or folder with images';
    if ~exist(rootDir,'dir')
        rootDir = pwd;
    end
    REFilter ='\.(tif(f?)|raw|prd)$'; %'\.tif(f?)$'
    outFile = uipickfiles('FilterSpec',rootDir,'REFilter',REFilter,'Prompt',prompt,'Output','char');
    outFile = strtrim(string(outFile));
    
    outFile = findAllRawFiles(outFile);
    [filePath,fileNames,ext] = fileparts(outFile);
    fileNames = fileNames+ext;
    filePath(isfolder(outFile))=outFile(isfolder(outFile));
    fileNames(isfolder(outFile)) = "";

    setpref(class(app),'lastPath',filePath(end));
        
    tempParam.pth = "";
    tempParam.imtoken = "";
    tempParam.wsfn = "";
    tempParam.nums = 1;
    tempParam.doRegistration = [];
    tempParam.BGRemovalRadius = [];
    % Process non-multiplane TIFFs
    paramsArray = repmat(tempParam, length(fileNames), 1);
    
    for i = 1:length(fileNames)
        paramsArray(i).imtoken = fileNames(i);
        paramsArray(i).pth = filePath(i);
        paramsArray(i).nums = 1;
    
        paramsArray(i) = updateRegAndBG(paramsArray(i));
        paramsArray(i) = updateWSFN(paramsArray(i));
    end
end


function param = updateRegAndBG(param)
    savefnTS = getFileName(param, 'TS');
    if exist(savefnTS,"file")
        fileContentsTS = load(savefnTS);
        if isfield(fileContentsTS,"didReg") && ~isempty(fileContentsTS.didReg)
            param.doRegistration = fileContentsTS.didReg;
        end
        if isfield(fileContentsTS,"BGRemovalRadius") && ~isempty(fileContentsTS.BGRemovalRadius)
            param.BGRemovalRadius = fileContentsTS.BGRemovalRadius;
        end
    else
        regInfoFn = getFileName(param, 'RegInfo');
        if exist(getFileName(param, 'RegBGMov'),"file") && exist(regInfoFn,"file")
            param.doRegistration = 1;
            regInfo = load(regInfoFn);
            if isfield(regInfo,'bgRadius')
                param.BGRemovalRadius = regInfo.bgRadius;
            end
        end
    end
end

function param = updateWSFN(param)

    matFile = dir(fullfile(param.pth, '*.mat'));

    if isscalar(matFile)
        param.wsfn = string(matFile.name);
        return
    end
    
    h5File = dir(fullfile(param.pth, '*.h5'));
    if isscalar(h5File)
        param.wsfn = string(h5File.name);
        return
    end
end

function files = findAllRawFiles(folderPaths)
    files = folderPaths(isfile(folderPaths));
    folders = folderPaths(~isfile(folderPaths));
    allFilesInFolders = arrayfun(@(a) dir(fullfile(a, '**', '*.raw')),folders,'UniformOutput',false);
    baseFolders = folders(cellfun(@isempty,allFilesInFolders));
    allFilesInFolders = cat(1,allFilesInFolders{:});
    if isempty(allFilesInFolders)
        allFilesInFolders = string.empty;
    else
        allFilesInFolders = string({allFilesInFolders.folder})'+repmat(filesep,[length(allFilesInFolders),1])+string({allFilesInFolders.name})';
    end
    files = cat(1,files,baseFolders,allFilesInFolders);
end