function [h,EventSpikes, ratevals, times, TrialReference,rastpos,histpos,spikecolor]= EventPethraPlot(spikes,events,varargin)
% a new and improved spike raster plotting function
%
%
%
%
% INPUTS:
%       spikes: vector of spike timestamps has to be column vector i think
%       events: a vector of events, can be concatenated multiple events,
%       etc. column vector as well
%       varargin:
%           'SecBefore': seconds before, either scalar or vector matching
%           the number of events, if you normalize, it will be whatever
%           that duration is as a percentage of the first epoch
%           'SecAfter': seconds after event, either scalar or vector
%           matching the number of events ** may want to cahnge to just 2
%           numbers
%           'EventID': vector matching events with 1:n different ID's.
%           Each id will (or will not be sorted and the raster will be a different color
%           and the line plot will be different.
%           'EventSort': 1, or 0 if you dont want to sort multiple events
%           'EventColor': RGB colors to match the event id's, if not ill
%           choose
%           'Epochs': either a column of stop times, or an N x 2 matrix
%           N events by 2 start/stop, for epochs you want to stretch, it
%           will be the same color as the spikes but just high gamma if
%           you dont separately designate spiketag
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
%           'AxisCoords': position in the graph you want your plot, this is
%           helpful if you want to add multiple of these to a single figure
%           'HistBins' make your own histogram bins, its a number of bins
%           so if you want it in seconds, figure out how many tot seconds
%           there are
%           'EventLines' A big horizontal line between your event rasters
%           'AddPatch' 1 or 0, 1 for a patch behind mean fr
% OUTPUTS:
%        h=figure handle
%        spikematrix: 3 col vector with trialidx, trialtime, and spiketag
%        lineplots: the x and y for the lineplots, will have multiple cells if
%        there are multiple epochs
% this functoin requires the following functions
% event_spikes, pe_th

% first get inputs in order

p=inputParser;

ValidateNumber= @(a) isnumeric(a);
ValidateTF= @(a) islogical(a);
addOptional(p,'SecBefore',1,ValidateNumber);
addOptional(p,'SecAfter',2,ValidateNumber);
addOptional(p,'EventID',1);
addOptional(p,'Epochs',0);
addOptional(p,'EventLines',0);
addOptional(p,'NormalizeTime',false);
addOptional(p,'Xscale',[]);
addOptional(p,'spiketag',[]);
addOptional(p,'EventSort',0,ValidateNumber);
addOptional(p,'AxisCoords',[]);
addOptional(p,'HistBins',30);
addOptional(p,'colormap',[]);
addOptional(p,'EventColors',[]);
addOptional(p,'AddPatch',0);
addOptional(p,'UseLabels',1);

addOptional(p,'NewMax',0,ValidateNumber);

parse(p,varargin{:});

SecBefore=p.Results.SecBefore;
SecAfter=p.Results.SecAfter;
EventID=p.Results.EventID;
Epochs=p.Results.Epochs;
Timecrunch=p.Results.NormalizeTime;
Xscaling=p.Results.Xscale;
spiketag=p.Results.spiketag;
mycolormap=p.Results.colormap;
EventColors=p.Results.EventColors;
EventSort=p.Results.EventSort;
AxisCoords=p.Results.AxisCoords;
histbins=p.Results.HistBins;
newmax=p.Results.NewMax;
addpatch=p.Results.AddPatch;
EventLines=p.Results.EventLines;
UseLabels=p.Results.UseLabels;
%% the algorithm;
% 1. get start, stop times, and their ID (what type of event they are)
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

% get stop times (also in seconds from 0)
if length(SecAfter)==1
    SecAfter=repmat(SecAfter,NumEvents,1);
elseif length(SecAfter)<NumEvents
    warning('too few SecAfter timestamps, filling in with mean time');
    SecAfter(end:NumEvents,1)=mean(SecAfter);
end

% epochs will be real timestamps
if any(Epochs)
    if size(Epochs,1)<NumEvents
        warning('epochs variable is too small, using set time after event start');
        % create our epochs
        Epochs(end:NumEvents,1)=mean(Epochs(:,1));
        Epochs(end:NumEvents,2)=mean(Epochs(:,2));
        if size(epochs,2)<2
            Epochs=[eventmat eventmat+(Epochs-eventmat(1))];
        else
            Epochs=[eventmat-(eventmat(1)-Epochs(1)) eventmat+(Epochs(1,2)-eventmat(1))];
        end
    end
    % now make sure we change our epochs to be relative to start times
    Epochs(:,1)=Epochs(:,1)-events; Epochs(:,2)=Epochs(:,2)-events;
else
    % if no epochs make it empty so we can add it to the list
    Epochs=[];
end

% SecBefore and SecAfter are both real now



%%%%%%%%% pull the spike times %%%%%%%%%%%%%

% for each event get cell of spike times, these will be used latah
for i=1:NumEvents
    % bracket spikes to be within our epoch
    [~,~,EventSpikes{i},inds{i}]=event_spikes(spikes,eventmat(i),SecBefore(i),SecAfter(i));
    % get these spikes to all be relative to zero
    %EventSpikes{i}=EventSpikes{i}-eventmat(i); updated event_spikes to fix
    %this
    % get our events to be relative to zero as well, were gonna have to
    % plot those out as well
    if any(spiketag) 
        spikecolor{i}=spiketag(inds{i}{1}); 
    else
        spikecolor=[];
    end
    evwindow(i,:)=[eventmat(i)-SecBefore(i) eventmat(i)+SecAfter(i)];
end

%%%%%%%%% see if we need to designate colors %%%%%%%%%%%%%%%%%

% how many event types do we have 1 is default, meaning one color.  if
% eventID starts with a zero we'll have
EventTypes=ceil(length(unique(EventID))/2)*2;
% fix colors to an even number, so we can increase contrast between consec
% events
if length(EventColors)<EventTypes && ~isempty(EventColors)
    fprintf('You didnt specify enough colors, needs to be row for each event type \n');
    EventColors=[];
end


    if length(EventID)<2
        allEventColors=repmat([0 0 0],NumEvents,1);
        EventID=ones(NumEvents,1);
    else
        
        % get our color palate
        % if we have our own, use it
        % sometimes not all colors are represented though
        if length(EventColors)>=EventTypes
            Colors=EventColors;
            % now attach to each trial
            allEventColors=Colors(EventID,:);
        % if we dont specify colors, or we dont have enough, make our own
        % colormap (this is a random set of colors, may need fixing)
        else
            % otherwise were gonna randomize it
            Colors=jet(2*ceil(EventTypes/2));
            %Colors=parula(2*ceil(EventTypes/2));
            
            % now assign each event a color
            allEventColors=zeros(NumEvents,3);
            % for each event type
            idxE=1; idxO=1;
            for i=1:EventTypes
                % if were at an even number, start half way through the rainbow
                
                if rem(i,2)==0
                    % so it goes 1,3,2,4 if there are 4 colors
                    % nah, we'll do just parula in good order
                    
                    for q=1:3, allEventColors(EventID==i,q)=Colors(idxE+(floor(EventTypes/2)),q); end
                    idxE=idxE+1;
                    % otherwise, start at beginning of rainbow
                else
                    
                    for q=1:3, allEventColors(EventID==i,q)=Colors(idxO,q); end
                    idxO=idxO+1;
                end
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
%fprintf('wait');
% lets define, or redefine EventID's;
% if we decided to order any events
if any(EventSort)
    % Start with unsorted events if we add 0 infront
    
    if EventSort(1)==0
        for i=1:length(EventID), unsorted(i)=~any(EventID(i)==EventSort); end
        % turn unsorted events into 0's, and add 1 to all events, that way
        % we can index our sort to replace the original sorting numbers
        if exist('unsorted','var'),EventID(unsorted)=0; EventID=EventID+1; end
        EventID=EventSort(EventID)';
    else
        % or put unsorted events at the end
        for i=1:length(EventID), unsorted(i)=~any(EventID(i)==EventSort); end
        if exist('unsorted','var'), EventID(unsorted)=max(EventSort+1); EventSort=[EventSort length(EventSort)+1]; end
        EventID=EventSort(EventID)';
    end
else
    % otherwise, reset the event ID's and make only one sort category
    EventSort=1; EventID=ones(NumEvents,1);
end


%%%%%%% Aggregate information for each trial %%%%%%%%%
% crucial info is start and stop for gray, epoch (which may be colored), and
% the colors of the tick marks, which epoch it belongs to

% 1 is event time (for ordering)
% 2,3 are start and end of trial,
% 4 is event type
% 5:7 are the color of the ticks/background
% 8:9 may be the epoch
TrialReference=[eventmat evwindow EventID allEventColors Epochs];
% sort first by trialtype, then by trial timestamp
if length(EventSort)>1
    [TrialReference,reorder]=sortrows(TrialReference,[4 1]);
    EventSpikes=EventSpikes(reorder);
    if iscell(spikecolor)
        spikecolor=spikecolor(reorder);
    end
end


%%%%%%%%% to plot it out %%%%%%%%%%%%%%%%%

%Prepare axes for plotting, the raster goes on top
if any(AxisCoords)
    h=get(gcf);
else
    h=figure;
    AxisCoords=get(gca,'Position');
end
% plot size ratio, 
plotratio=.3; % what percent do you want the histogram to be
rastpos=subplot('position',[AxisCoords(1) AxisCoords(2)+AxisCoords(4)*plotratio ...
    AxisCoords(3) AxisCoords(4)*(1-plotratio)]);
histpos=subplot('position',[AxisCoords(1) AxisCoords(2) ...
    AxisCoords(3) AxisCoords(4)*plotratio]);


% this is a place holder to colorcode spikes by their phase
spikereference=[];
if iscell(spikecolor)
    %get all data and make our bounds based on that
    allphases=cell2mat(spikecolor');
    if any(allphases) % if there are spikes to sort
        bounds=linspace(min(allphases),max(allphases),128);
        cmap=hsv(128);
        % now for each trial slot the color of the spike into that
        for i=1:length(spikecolor)
            [~,binmap]=histc(spikecolor{i},bounds);
            try % if there are nans, we gotta nan these guys out
                spikereference{i}=cmap(binmap,:);
            catch
                spikereference{i}=ones(length(binmap),3);
            end
        end
    end
end
axeslabels=UseLabels;

MakeRaster(rastpos,TrialReference,EventSpikes,spikereference,EventLines);
set(rastpos,'ydir','reverse',...
    'xlim',[-max(SecBefore) max(SecAfter)],...
    'ylim',[1 NumEvents+1]);
set(rastpos,'box','on');
if axeslabels, ylabel('Trial'); end

[ratevals,times]=MakeHisto(histpos,TrialReference,spikes,[SecBefore(1) SecAfter(1)],histbins,newmax,addpatch);
box on;

if axeslabels
    xlabel('Seconds From Sample Onset');
    ylabel('Rate, Hz');
end

end

function [ratevals,times]=MakeHisto(histhandle,TrialReference,CellSpikes,window,histbins,newmax,addpatch)
% address the plot were gonna use
subplot(histhandle); hold on
% for each event type, get after the trial types
TrialIDs=unique(TrialReference(:,4));
for i=1:length(TrialIDs)
    clear mytrials;
    % get the trials that you want to pull (for each color)
    mytrials=TrialReference(TrialReference(:,4)==TrialIDs(i),:);
    % now run peth to find where the spikes are
    [xtimes(i,:),rates(i,:),mymax(i)]=pe_th(CellSpikes, mytrials(:,1), window(1), window(2),histbins);
    plot(xtimes(i,:),(rates(i,:)),'color',mytrials(1,5:7),'LineWidth',3);
    hold on
    
    if addpatch
         binsize=(window(1)+window(2))/histbins;
         mybins=-window(1)+binsize/2:binsize:window(2);
        try
           errmat=[];
            for j=1:length(mybins)
                [~,errmat(:,j)]=event_spikes(CellSpikes,mytrials(:,1)+mybins(j),-binsize/2,binsize);
            end
            errs=SEM(errmat,1);
            mypatch=patch([ xtimes(i,:) fliplr(xtimes(i,:))],...
                [rates(i,:)+errs fliplr(rates(i,:))-fliplr(errs)],mytrials(1,5:7));
            set(mypatch,'EdgeAlpha',0,'FaceAlpha',.2);
        catch
        end
    end
            
end
ratevals=rates;
times=xtimes;

% if no spikes, kill figure
if max(mymax)==0
    if newmax>0
        set(gca,'ylim',[0 max(newmax)]);
    end
    %close gcf
    return
    % otherwise make some lines, or make a patch
else
    if newmax>0
        mymax=max([newmax mymax]);
    end
    set(gca,'xlim',[-window(1) window(2)],...
        'ylim',[0 max(mymax)*1.2]);
    % if length(TrialReference,2)>7
    
    shadow=[0 8];
    % else
    if ~isempty(shadow)
    line([shadow(1) shadow(1)],[0 100],'Color','r');
    line([shadow(2) shadow(2)],[0 100],'Color','r');
    % or patch (use trialreference)
    end
end
end

function MakeRaster(rasterhandle, TrialReference,EventSpikes,spikereference,EventLines)
% eventspikes is an n cell aray where n is trials, in each cell is the ts
% for each spike (as it relates to the event
% trialreference is either a 3 x n vector, n is the trial, or its an n cell
% array where each cell is a trial. in each cell its a 3 x n spikes per
% trial
%%%%%%%%%% SET LINE WIDTH HERE %%%%%%%%%%%%%
linewidth=.2;

% will have already made a position for this plot
subplot(rasterhandle);
set(gca,'Xtick',[]);
% first patch the whole raster if youre coloring your spikes
if ~isempty(spikereference)
    xpatch=[TrialReference(1,2:3) flip(TrialReference(1,2:3),2)]-TrialReference(1);
    ypatch=[0 0 length(TrialReference)+1 length(TrialReference)+1];
    patch(xpatch*1.1,ypatch,[.95 .95 .95],'LineStyle','none');
end
% now fill in each trial


for k=1:length(TrialReference(:,1))
    
    % now put in the spikes and colorcode them
    ThisRaster=EventSpikes{k};
    if any(ThisRaster)
        for q=1:length(ThisRaster)
            linehandle = line([ThisRaster(q),ThisRaster(q)],[k k+.95]);
            if ~isempty(spikereference)
                set(linehandle,'Color',spikereference{k}(q,1:3),'LineWidth',linewidth);
            else
                set(linehandle,'Color',TrialReference(k,5:7),'LineWidth',linewidth);
            end
        end
    end
    % first gray in the whole trial out
    if isempty(spikereference)
        xpatch=[TrialReference(k, 2:3) flip(TrialReference(k,2:3),2)]-TrialReference(k,1);
        ypatch=[k k k+1 k+1];
        thispatch=patch(xpatch, ypatch, [.95 .95 .95],'LineStyle','none');
        alpha(thispatch,.7);
    end
    
    % if theres an epoch, put that in too
    if size(TrialReference,2)>7
        xpatch=[TrialReference(k, 8:9)  flip(TrialReference(k,8:9))];
        ypatch=[k k k+1 k+1];
        % color it the right color
        %thispatch=patch(xpatch, ypatch, TrialReference(k,5:7),'LineStyle','none');
        % or just gray
        thispatch=patch(xpatch, ypatch, [1 .7 .7],'LineStyle','none');
        alpha(thispatch,.1);
    end  
end
% put a pair of lines to designate the delay
shadow=[0 8];
if ~isempty(shadow)
    line([shadow(1) shadow(1)],[0 length(TrialReference(:,1))+1],'Color','r');
    line([shadow(2) shadow(2)],[0 length(TrialReference(:,1))+1],'Color','r');
end
if EventLines==1
    [a]=find(diff(TrialReference(:,4))~=0);
    xvals=TrialReference(k, 2:3)-TrialReference(k,1);
    for k=1:length(a)
        line(xvals,[a+1 a+1],'Color','k');
    end
end


%if ~isempty(spikereference); set(rasterhandle,'Color',[.2 .2 .2]); hold on; end
end







