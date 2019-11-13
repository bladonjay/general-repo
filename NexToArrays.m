function [unitTS,eventTS,units,events,filename]=NexToArrays(filename, dontask)
%This script reads a *.NEX file and outputs an arrays of timestamps for 
% each identified unit and event.

%Input: filename = name of *.nex file to be analyzed; if no variable is
%   entered a GUI will appear
%       dontask = binary variable enabling query to extract data from file when
%       dontask == 1  i.e. if you enter a var thats not 1 for dontask 
%       you will be asked whether this is the right file
%Output: unitTS = Unit timestamps in an array of matrices
%       eventTS = Event timestamps in an array of matrices
%       units = Unit names
%       events = Event names

if nargin ~= 2
    dontask=1;
end
if nargin < 1
    [filename,pathname]=uigetfile('*.nex', 'Select a Nex file');
    filename=[pathname filename];
end

fclose('all');
[nvar, names, types] = nex_info(filename);
name2={};
for m=1:nvar
    name2=[name2; names(m,:)];
end
units = name2(types==0);
events = name2(types==1);
disp(['There are ', num2str(length(units)), ' unit(s) and ', num2str(length(events)), ' event type(s) in this file.'])
if dontask ~= 1
    response=input('Would you like to extract timestamp data from this file? [Y/N]  ', 's');
    if ~strcmpi(response,'Y')
        return
    end
end
unitTS={};
eventTS={};

% changed to display each event name.  nex_ts displays the filename for
% each output
% there are never empty units, so we can run all units
for m=1:length(units)
    fprintf(' unit %d \n', m);
    [n, ts] = nex_ts(filename, units{m});
    unitTS{m}=ts;
end
keep=[];
for m=1:length(events)
    disp(events{m});
    [n, ts] = nex_ts(filename, events{m});
    eventTS{m}=ts;
    keep(m)=~isempty(ts);
end

% now delete all the empty events and units
events=events(logical(keep));
eventTS=eventTS(logical(keep));
end
