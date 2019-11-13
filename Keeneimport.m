function [unitTS,eventTS,units,events,filename]=Keeneimport(dirname)

if nargin~=1
    dirname=uigetdir;
end


if ispc
    filenames=dir([dirname '\*.nex']);
    % backslash for pc file systems
    filenames=struct2cell(filenames);
    filenames=filenames(1,:);
    nexfiles=strcat([dirname '\'],filenames);
end
if ismac
    filenames=dir([dirname '/*.nex']);
    % backslash for pc file systems
    filenames=struct2cell(filenames);
    filenames=filenames(1,:);
    nexfiles=strcat([dirname '/'],filenames);
end

unitTS={};
eventTS={};
units={};
events={};
filename={};

for m=1:numel(nexfiles)
    [curunitTS,cureventTS,curunits,curevents,curfilename]=NexToArrays(nexfiles{m}, true);
    curunitTSa=cell(length(curunitTS),1);
    for n=1:length(curunitTS)
        curunitTSa{n}=unique(curunitTS{n}); %This gets rid of repeat TS.
        if length(curunitTSa{n})~=length(curunitTS{n})
            warning('Repeat spike timestamps. Deleting repeats.');
        end
    end
    unitTS=[unitTS; curunitTSa];
    eventTS=[eventTS; cureventTS'];
    units=[units; curunits];
    events=[events; curevents];
    filename=[filename; curfilename];
end

end
