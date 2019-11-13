function [unitTSsub,eventTSsub,units,events,filename]=NexToMultiArrays(starts,stops)
[unitTS,eventTS,units,events,filename]=NexToArrays;

%'starts' and 'stops' are 1D arrays of start and stop timestamps dividing a
%*.nex file


%Had to find 'start' and 'stop' events circuitously because strcmp didn't
%like the Nex format.
Xstartind=strfind(events,'Start');
startind=[];
for m=1:length(Xstartind)
    startind=[startind;isempty(Xstartind{m})];
end
startind=find(startind==0,1);

Xstopind=strfind(events,'Stop');
stopind=[];
for m=1:length(Xstartind)
    stopind=[stopind;isempty(Xstopind{m})];
end
stopind=find(stopind==0,1);
    
if isempty(starts) || isempty(stops)
    if length(eventTS{startind})~=length(eventTS{stopind})
        disp('Unequal Number of Start and Stop flags');
        return
    end
    starts=eventTS{startind};
    stops=eventTS{stopind};
end

%Preallocate
unitTSsub=cell(length(starts),1);
for m=1:length(unitTSsub)
    unitTSsub{m}=cell(length(unitTS),1);
end

eventTSsub=cell(length(starts),1);
for m=1:length(eventTS{startind})
    eventTSsub{m}=cell(length(eventTS)-2,1);
end

%Split Units in Nex-file
for m=1:length(unitTSsub)
    for n=1:length(unitTS)
        curUnitTS=unitTS{n}(unitTS{n}>=starts(m) & unitTS{n}<stops(m));
        curUnitTS=curUnitTS-starts(m);
        unitTSsub{m}{n}=curUnitTS;
    end
end

%Split Events in Nex-file
XeventTS=eventTS;
XeventTS([startind stopind])=[];
events([startind stopind])=[];

for m=1:length(eventTSsub)
    for n=1:(length(eventTS)-2)
        curEventTS=XeventTS{n}(XeventTS{n}>=starts(m) & XeventTS{n}<stops(m));
        curEventTS=curEventTS-starts(m);
        eventTSsub{m}{n}=curEventTS;
    end
end