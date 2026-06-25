function analysisDir = getAnalysisDir(params)

analysisDir = fullfile(params.pth, 'Analysis');
if ~exist(analysisDir, 'Dir')
    mkdir(analysisDir);
    disp(['Created directory ' analysisDir '!!!']);
end