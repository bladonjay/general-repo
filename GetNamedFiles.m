function [fileList] = GetNamedFiles(dirName,stringmatch,inputdlg)
% pull m files with specific name. stringmatch is obligatory, dirName and
% inputdlg arent

% JH BLaodn


if ~exist('inputdlg','var')
    inputdlg='choose directory with your files!';
end
if ~exist('dirName','var')
    dirName=uigetdir('',inputdlg);
elseif isempty(dirName)
    dirName=uigetdir('',inputdlg);
end
% pull the spksevs files from the folder
fileList = getAllFiles(dirName);
session=cellfun(@(a)  any(regexpi(a,stringmatch))   ,fileList,'UniformOutput',0);
session=cell2mat(session);
fileList=fileList(session);

end

