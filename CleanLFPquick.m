function [LFPclean, LFPts] = CleanLFPquick(LFPstruct)
% Function CleanLFPquick stitches together the LFP, assigns timestamps to
% each datapoint, and removes the datapoints in which the signal has
% railed.

% In the future, I ought to filter out low frequency noise



% first lets reconstitute our LFP;
fragments=LFPstruct.fragmentStarts;
timefrags=LFPstruct.timestamps;
% add first timestamp
LFPstruct.time=[];
dt=1/LFPstruct.ADFrequency;
% last fragment we have to treat different because it doesnt have a stop
% time
for i=1:length(fragments)-1
    % find number of datapoints recorded
    thislength=fragments(i+1)-fragments(i);
    % now step forward till right before the next frag starts
    fragts=(dt:dt: thislength *dt) + timefrags(i);
    LFPstruct.time=[LFPstruct.time fragts];
end


% adress last fragment, unique because now we account for last timestamp:
lastlength=length(LFPstruct.data)-fragments(end)+1;
if length(fragments)>1
    lastfrag=(dt:dt:(lastlength*dt)) + fragts(end);
else
    lastfrag=(dt:dt:(lastlength*dt));
end  
LFPstruct.time=[LFPstruct.time lastfrag];

ts=LFPstruct.time'; data=LFPstruct.data;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% filter the LFP %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scrub timestamps and data that are railing
% I'll use the max and min, so either i will scrub two datapoints, or all
% those that rail

% So then the data is all concatenated, and where there are borders, i can
% go back with my timestamps and remove

bottomrail=min(data);
removeind=find(data==bottomrail);
data(removeind)=[]; ts(removeind)=[];

toprail=max(data);
removeind=find(data==toprail);
data(removeind)=[]; ts(removeind)=[];

LFPclean=data; LFPts=ts;
% this should yield a cleaned data trace
end

