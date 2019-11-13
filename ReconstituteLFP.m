function [LFP] = ReconstituteLFP(contvars)
% [LFP]=ReconstituteLFP(contvars)
% LFP data come in as fragments because they have a lot of ts, so you have
% to rebuild the timestamps from the datset.  I would advise not doing this
% for every LFP because your data file will be HUGE.

if iscell(contvars)
    for i=1:length(contvars)
        LFP(i).name=contvars{i}.name;
        LFPstruct=LFPfrag2cont(contvars{i});
        LFP(i).ts=LFPstruct.time;
        LFP(i).data=LFPstruct.data;
    end
elseif isstruct(contvars)
    for i=1:length(contvars)
        LFP(i).name=contvars(i).name;
        LFPstruct=LFPfrag2cont(contvars(i));
        LFP(i).ts=LFPstruct.time;
        LFP(i).data=LFPstruct.data;
    end
end
end



function [LFPstruct]=LFPfrag2cont(contvar)
fragmentStarts=contvar.fragmentStarts;
timestamps=contvar.timestamps;
dt=1/contvar.ADFrequency;


% add first timestamp
LFPstruct.time=[];
LFPstruct.data=contvar.data';
fragts=[];
% last fragment we have to treat different because it doesnt have a stop
% time
for i=1:length(fragmentStarts)-1
    % find number of datapoints recorded
    thislength=fragmentStarts(i+1)-fragmentStarts(i);
    % now step forward till right before the next frag starts
    fragts=(dt:dt: thislength *dt) + timestamps(i);
    LFPstruct.time=[LFPstruct.time fragts];
end

% adress last fragment, unique because now we account for last timestamp:
lastlength=length(LFPstruct.data)-fragmentStarts(end)+1;
lastfrag=(dt:dt:(lastlength*dt)) + timestamps(end);
LFPstruct.time=[LFPstruct.time lastfrag];

end