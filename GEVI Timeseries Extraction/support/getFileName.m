function fn = getFileName(params, filetype)

analsisDir = getAnalysisDir(params);
switch filetype
    case 'ExpInfo'
        fn = 'ExpInfo(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'WSData'
        fn = 'WSData(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'ROI'
        fn = 'ROI(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'regMask'
        fn = 'regMask(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'RegInfo'
        fn = 'RegInfo(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'Mov'
        fn = 'Mov(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'RegBGMov'
        fn = 'RegBGMov(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.tif';
    case 'StimDFFMov'
        fn = 'StimDFFMov(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.tif';
    case 'Movavg'
        fn = 'Movavg(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.tif';
    case 'corrim'
        fn = 'corrim(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.tif';
    case 'dfRangeIm'
        fn = 'dfRangeIm(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.tif';
    case 'dfRatioim'
        fn = 'dfRatioim(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.tif';
    case 'Disp'
        fn = 'Disp(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'TS'
        fn = 'TS(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'Snip'
        fn = 'Snip(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'OptoData'
        fn = 'OptoData(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'BootstrapInfo'
        fn = 'BootstrapInfo(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'GLM'
        fn = 'GLM(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
    case 'DMDPatternInfo'
        fn = 'DMDPatternInfo(' + params.imtoken + ')_' + num2str(params.nums(1)) + '-' + num2str(params.nums(end)) + '.mat';
        if exist(fullfile(analsisDir,fn),'file')~=2
            fn = 'DMDPatternInfo.mat';
        end
        if exist(fullfile(analsisDir,fn),'file')~=2
            fn = fullfile('..','..','Analysis','DMDPatternInfo.mat');
        end
    case 'wsfn'
        fn = "";
        if isfield(params,"wsfn") && (exist(params.wsfn,"file") ==2 || exist(fullfile(params.pth,params.wsfn),"file") ==2)
            if exist(params.wsfn,"file") ==2 %2 is file (7 is folder)
                fn = params.wsfn;
            elseif exist(fullfile(params.pth,params.wsfn),"file") ==2  %2 is file (7 is folder)
                fn = fullfile(params.pth,params.wsfn);
            end
        end
        return
        
    otherwise
        throw(['Could not find filename for ' filetype '!!!!']);
end

fn = fullfile(analsisDir,fn);

if ~exist(fn,'file') && exist("fn2",'var') && exist(fn2,'file')
    fn = fn2;
end

