function [unitTS,unitnames,filename]=Keeneimportunit();
if nargin~=1
    dirname=uigetdir;
end

filenames=dir([dirname '\*.nex']);
filenames=struct2cell(filenames);
filenames=filenames(1,:);
nexfiles=strcat([dirname '\'],filenames);

unitTS={};
unitnames={};
filename={};

for m=1:numel(nexfiles)
    [curunitTS,cureventTS,curunitnames,curevents,curfilename]=NexToArrays(nexfiles{m}, true);
    curunitTSa=cell(length(curunitTS),1);
    for n=1:length(curunitTS)
        curunitTSa{n}=unique(curunitTS{n}); %This gets rid of repeat TS.
        if length(curunitTSa{n})~=length(curunitTS{n})
            warning('Repeat spike timestamps. Deleting repeats.');
        end
    end
    unitTS=[unitTS; curunitTSa];
    unitnames=[unitnames; curunitnames];
    filename=[filename; curfilename];
end