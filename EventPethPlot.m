[h,spikematrix, lineplots]= EventPethPlot(spikes,events,varargin);
% a new and improved spike raster plotting function
%
%
%
%
% INPUTS:
%       spikes: vector if spike timestamps
%       events: a vector of events, can be concatenated multiple events,
%       etc.
%       varargin:
%           'SecBefore': seconds before, either scalar or vector matching
%           the number of events, if you normalize, it will be whatever
%           that duration is as a percentage of the first epoch
%           'SecAfter': seconds after event, either scalar or vector
%           matching the number of events
%           'EventID': vector matching events with 1:n different ID's.
%           Each id will (or will not be sorted and the raster will be a different color
%           and the line plot will be different.
%           'EventSort': 1, or 0 if you dont want to sort multiple events
%           'EventColor': RGB colors to match the event id's, if not ill
%           choose
%           'Epochs': an N x 2 matrix N events by 2 start/stop,
%           for epochs you want to stretch, it will be the same color as
%           the spikes but just high gamma if you dont separately designate
%           spiketag
%           'NormalizeTime': if you want to normalize time, time will be
%           from 0(SecBefore) to 1(SecAfter) and everything will be
%           squished to be inside that (or start epoch to stop epoch plus
%           and minus your secbefore and secafter
%           'spiketag': a tag for the color of each spike- will overwrite
%           the color of your spikes if youre using eventid.  This would be
%           mostly for the phase of a spike to the lfp or its relationship
%           to another cell maybe
%           'colormap': a string, it could be jet, hsv, parula, or
%           whatever, the standard is hsv
% OUTPUTS:
%        h=figure handle
%        spikematrix: 3 col vector with trialidx, trialtime, and spiketag
%        lineplots: the x and y for the lineplots, will have multiple cells if
%        there are multiple epochs


% first get inputs in order
p=inputParser;

ValidateNumber= @(a) isnumeric(a);
ValidateTF= @(a) islogical(a);
addOptional(p,'SecBefore',0,ValidateNumber);
addOptional(p,'SecAfter',1,ValidateNumber);
addOptional(p,'EventID',1);
addOptional(p,'Epochs',0);
addOptional(p,'NormalizeTime',false);
addOptional(p,'Xscale',[]);
addOptional(p,'spiketag',[]);
addOptional(p,'EventSort',0,ValidateNumber)
addOptional(p,'colormap');

parse(p,varargin{:});

SecBefore=p.Results.SecBefore;
SecAfter=p.Results.SecAfter;
EventID=p.Results.EventId;
Epochs=p.Results.Epochs;
Timecrunch=p.Results.NormalizeTime;
Xscaling=p.Results.Xscale;
spiketags=p.Results.spiketag;
mycolormap=p.Results.colormap;
EventSort=p.Results.EventSort;


%% the algorithm;
% 1. get start, stop times, and their ID
%       ~ if you have an epoch, secafter is seconds after your epochs 
% 2. pull spike data in for each trial
%       ~ for each spike data cell, you'll have a color index as well
% 3. pull epoch data in for each trial
%       ~ and designate your color scheme
% 4. crunch if you need to
%       ~ crunch by smapping each point as a percentage of the average
%       time, or part per thousand, not quite sure yet
% 5. plot it all out in raster and line plot below
%       ~ if theres a epoch, put that in first with a high alpah
%       ~ then plot each timestamp with color scheme, if you decide not to
%       reorder, you'll have to recombine your trials with their
%       corresponding color indices



%%%%%% first get the length of your whole epoch %%%%%%%%%%%

% first start times
eventmat=events(:,1);
NumEvents=length(eventmat);

% get start times (in seconds from 0)
if length(SecBefore)==1
    SecBefore=repmat(SecBefore,NumEvents,1);
elseif length(SecBefore)<NumEvents
    warning('too few SecBefore timestamps, filling in with mean time');
    SecBefore(end:NumEvents,1)=mean(SecBefore);
end

% get stop times (also in seconds from event)
if length(SecAfter)==1
    SecAfter=repmat(SecAfter,NumEvents,1);
elseif length(SecBefore)<NumEvents
    warning('too few SecAfter timestamps, filling in with mean time');
    SecAfter(end:NumEvents,1)=mean(SecAfter);
end

% epochs will be real timestamps
if any(epochs)
    if length(epochs,2)<NumEvents
        warning('epochs variable is a scalar, using set time after event start');
        Epochs=[eventmat-(eventmat(1)-Epochs(1)) eventmat+(Epochs(1,2)-eventmat(1))];
    end 
end
% now turn them into real timestamps so we can pull spike timestamps
if any(Epochs);
        SecBefore=(eventmat-Epochs(:,1))+SecBefore;
        SecAfter=Epochs(:,2)+SecAfter;
        % the kicker, this is fed into the event_spikes function
    
else
    SecBefore=eventmat-SecBefore;
    SecAfter=eventmat+SecAfter;
    % again, kicker, feed this into event_spikes fcn
    
end

%%%%%%%%% pull the spike times

% for each event get cell of spike times, these will be used latah
for i=1:NumEvents
    % bracket spikes to be within our epoch
    [evspikes{i}]=event_spikes(spikes,eventmat(i),SecBefore(i),SecAfter(i));
    % get these spikes to all be relative to zero
    evspikes{i}=evspikes{i}-eventmat(i);
    evwindow{i}=[SecBeforeReal
end

%%%%%%%%% see if we need to designate colors

EventTypes=length(unique(EventID));
% fix colors to an even number, so we can increase contrast between consec
% events
if length(EventID)<2
    EventColors=repmat([1 1 1],NumEvents,1);
    EventID=ones(NumEvents,1);
else
    % first create a rainbow
    Colors=hsv(2*ceil(EventTypes/2));
    % now assign each event a color
    EventColors=zeros(3,NumEvents);
    for i=1:length(EventID)/2
        % if were at an even number, start half way through the rainbow
        if rem(i,2)==0, EventColors(EventID==i,:)=Colors(i*2);
        % otherwise, start at beginning of rainbow 
        else EventColors(EventID==i,:)=Colors(i);
        end
    end
end


%%%%%%%%%%% CRUNCH, NORMALIZE time here %%%%%%%%%%%%%%
%%% i'll do that later


%%%%%%%%% sort events if were doing that %%%%%%%%%

% basically any number you add to sortevents will isolate that event, and
% it will do it in order, e.g. if you decide to sort event type 1 first,
% then its 1.  if you decide to sort event type 1 out, but put it at the
% end, type 0, 1

% lets define, or redefine EventID's;
% if we decided to order any events
if any(EventSort)
    % turn all unsorted events to last if we add 0 infront
    if EventSort(1)==0
        for i=1:length(EventID), unsorted(i)=~any(EventID(i)==EventSort); end
        % turn unsorted events into 0's
        EventID(unsorted)=0;
    else
        % or turn unsorted events into a new number
        for i=1:length(EventID), unsorted(i)=~any(EventID(i)==EventSort); end
        EventID(unsorted)=max(EventSort+1);
    end
end




%%%%%%%%% to plot it out %%%%%%%%%%%%%%%%%

%Prepare axes for plotting
if any(strcmpi('axeshandle',namevars))
    h.fig=get(axeshandle,'Parent');
    h.axes=axeshandle;
else
    h.fig=figure;
    h.axes=axes;
end
% for each trial type
for i=1:length(EventSort)
    TheseRasters=evspikes{EventID==EventSort(i)};
    TheseColors=EventColors(EventID==EventSort(i));
    % ignore epochs for now
    TheseEpochs=[];
    
    
end



    

