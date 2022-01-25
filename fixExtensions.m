function [newname,oldname] = fixExtensions(filename,varargin)
% function [outputArg1,outputArg2] = fixExtensions(filename,varargin)
% function fixExtensions allows you to fix the file extensions in a file
% it by default replaces all periods with underscores except for the last one, but you can
% choose to use the first if you want
% INPUTS: filename: full file path of file
%   OPTIONAL:
%        'Ext', number: number order of extension you want to use.  default
%        is 0, which is the last extension in the file

p = inputParser;
addOptional(p,'Ext',0); 
parse(p,varargin{:});
Ext=p.Results.Ext;
%[fname,fdir]=uigetfile();
[folderstring,namestring,ext]=fileparts(filename);

oldname=fullfile(folderstring,[namestring ext]);
namedots=find(filename(find(filename==filesep,1,'last')+1:end)=='.');
if length(namedots)>1
    switch Ext
        case 0
             namestring(namestring=='.')='_';
             newname=fullfile(folderstring,[namestring ext]);
        case 1
            newname=fullfile(folderstring,namestring(1:namedots(2)-1));
        case inf
            newname=fullfile(folderstring,[namestring(1:namedots(1)-1) ext]);

    end
    movefile(oldname,newname);
else
    newname=oldname;
end



end

