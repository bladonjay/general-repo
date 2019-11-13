function [data, nex] = jayNEX2CMBbatch(dirName,varargin)
%This is a simple wrapper function to batch process a folder containing various
%nex files into a CMBHome format.
% [CMBobject, nexstruct]= jayNEX2CMBbatch(directoryname)
% directory can contain one or many nex files, the nex files can contain
% flags, units, or both.  The output will create an aggregate nex struct as
% well as an aggregate CMBobject.  
% CALLS: MakeNEXMaster and NEX2CMBjay
% JHB 6-10-14
if ~exist('varargin','var') || isempty(varargin)
    varargin=0;
else
    varargin=varargin{1};
end
% first I read all the nex files and create an aggregate struct
if ~exist('dirName','var') || isempty(dirName)
    dirName=uigetdir('','choose a folder of nex files');
end

[nex,source]=MakeNexMaster(dirName);
%%%%%%%
%%to add the xy data
%%%%%%%
files=dir(dirName);
fileList={files.name};
fileList=fileList(cellfun(@any,strfind(fileList,'.mat')));
for i=1:length(fileList)
    % make sure we dont try to load the CMB object, because that sucks
    if isempty(strfind(fileList{i},'Session')) && isempty(strfind(fileList{i},'SpksEvs'))
        if ispc
            load([dirName '\' fileList{i}]);
        elseif ismac
            load([dirName '/' fileList{i}]);
        end
    end
end
nex.coords=cleanxydata;

[data, nex]= NEX2CMBjay(nex,varargin);

if ispc
    data.name= [source(max(strfind(source,'\'))+1:end) 'CMB.mat'];
    
    saver=input('would you like to save the file? 1 if yes 0 if no  ');
    if saver==1
        save([dirName '\' data.name(1:end-3) 'CMB.mat'],'data','nex');
    end
end

if ismac
    data.name= [source(max(strfind(source,'/'))+1:end) 'CMB.mat'];
    
    saver=input('would you like to save the file? 1 if yes 0 if no  ');
    if saver==1
        save([dirName '/' data.name(1:end-3) 'CMB.mat'],'data','nex');
    end
end

end

