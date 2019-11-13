%CleanWorkSpace
% CleanWorkSpace(varargin)


varlist=whos;
% grab a list of varnames
for i=1:length(varlist), varnames{i}=varlist(i).name; end


% if were leaving only a few vars
if exist('LeaveOut','var')
    % remove
    for i=1:length(varlist), killvars(i)=~any(contains(varnames{i},LeaveOut)); end
    % or if were removing only a few vars
elseif exist('LeaveIn','var')
    % otherwise, list em all out
    for i=1:length(varlist), killvars(i)=any(contains(varnames{i},LeaveIn)); end
else
    if length(varlist)>40
        varchunks=round(linspace(1,length(varnames),ceil(length(varnames)/40)+1));
        killvars=[];
        for i=1:length(varchunks)-1
            killvars=[killvars checkBox(varnames(varchunks(i):varchunks(i+1)))+varchunks(i)-1];
        end
    else
        killvars=checkBox(varnames);
    end
end
clear(varnames{killvars})