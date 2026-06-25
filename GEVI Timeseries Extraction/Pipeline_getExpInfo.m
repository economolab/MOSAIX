function Pipeline_getExpInfo(params)
fn.pth = params.pth;

savefn = getFileName(params, 'ExpInfo');
save(savefn, 'fn'); 

disp(' ');
disp('Saved experiment file: ');
disp(['filename: ' savefn]);















