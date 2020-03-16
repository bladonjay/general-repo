function [ratemap,topplot,finalcolormap,noccup2,cellspikes,velcoords] = SmoothPlacePlot2(coorddata,unit,varargin)
%FUNCTION [ratemap, rawratemap]= cell_SmoothPlacePlot(session, unit, varargin)
% This function creates a smoothed place plot of cell u.  There are two
% subplots, on top a trajectory overlaid with spike locations, and on
% bottom there is a heat plot showing the firing fields.  This can be
% combined with the 'epoch session' function if there are multiple epochs
% in the session.

% INPUTS;
% session is a session struct with a coords or edit_coords field
% unit is a struct of a cell: has to have a TS field with timestamps
% varargin is below, these are parameters to play with the plots.  the most
% salient are gaussdev (the standard deviation of the gaussian kernel) and
% the factor (the number of pixels per bin)
% the remaining variables are mostly constraints on data
%
% OUTPUTS
% This function will create a figure, with the handle g
% ratemap is a grid matrix where each cell is the firing rate at that pixel
% rawratemap is the ratemap with no smoothing

% notes;
% this smoothes after spikes are divided, like the hasselmo lab, unlike the
% mosers lab
% The gaussian kernel size can be as large as you want
% There is a variable at line
% this is also pretty messy, I should fix this in the future

%
% JHB 9-5-14
% JHB 3-9-20

p=inputParser;
addOptional(p,'Gausswin',20); % in pixels of final map
addOptional(p,'Gaussdev',1); % in pixels of final map
addOptional(p,'PixPerCM',6); % pixels per bin
addOptional(p,'CMPerBin',2); % bin as in pixel for rate map
addOptional(p,'minvelocity',5); % cm/sec
addOptional(p,'mintimespent',.2);  % or about 5 timestamps
addOptional(p,'minvisits',2); % visits to that bin
addOptional(p,'velsmooth',.2); % in seconds
addOptional(p,'spikethreshold',50); % nspikes you have to have to do the expt
addOptional(p,'clim',[]); 
addOptional(p,'EpochName','Whole Session');
addOptional(p,'RealLim',1); 
addOptional(p,'Grayout',0);
addOptional(p,'Timejumpdelay',.3); % how many seconds is the max to comb over
addOptional(p,'suppress',0);
addOptional(p,'onlyheat',0);
addOptional(p,'realcoords',0);
addOptional(p,'scalecoorddata',false);
addOptional(p,'VelocityData',[]);
parse(p,varargin{:});

% if we already have smoothed velocity data
VelData=p.Results.VelocityData;

%%%%% these are plotting parameters %%%%%%%
spikethreshold=p.Results.spikethreshold; % fewest spikes the cell has to have to plot
Gausswin=p.Results.Gausswin; % size of the gaussian kernel
Gaussdev=p.Results.Gaussdev; % STD of gaussian kernel
PixPerCM=p.Results.PixPerCM; % How many raw pixels per cm
CMPerBin=p.Results.CMPerBin; % how many cm per bin for final colormap?
clim=p.Results.clim; % if you want to scale the rateplot based on a set firing rate max, e.g.
RealLim=p.Results.RealLim; % use max firing rate or max reasonable?
grayout=p.Results.Grayout; % gray out no spikes?
EpochName=p.Results.EpochName; % Name the epoch for the figure
suppress=p.Results.suppress; % suppress plot
onlyheat=p.Results.onlyheat; % only heat plot
scalecoords=p.Results.scalecoorddata; % adjust coords to fit video?


%%%%% these are cutoff parameters %%%%%%%%
minvelocity=p.Results.minvelocity; % velocity cutoff to remove any stationary time
mintimespent=p.Results.mintimespent; % min time spent in each pixel
minvisits=p.Results.minvisits; % min visits to the pixel
velsmooth=p.Results.velsmooth; % time in seconds to smooth the velocity
timejumpdelay=p.Results.Timejumpdelay; % time delay in seconds that is too long to interpolate
realcoords=p.Results.realcoords; % real coordinates so that you can compare pixels across maps


% output variables;
ratemap=[];
rawratemap=[];
topplot=[];
finalcolormap=[];
noccup2=[];
cellspikes=[];
velcoords=[];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% make sure ts always goes up and no repeats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,uniquets]=unique(coorddata(:,1));
coorddata=coorddata(uniquets,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get velocity and make velocity filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if our coords arent the same as the video file
if scalecoords
    coorddata(:,2:end)=coorddata(:,2:end)*0.6246;
end

dt = diff(coorddata(:,1)); % difference between consecutive coords
% chop off the first timestamp because we dont know its elapsed time
timestamps = coorddata((2:end),1); % remove coords so dt1 is for before tcoord1
xycoords = coorddata((2:end),2:3); % save raw xy data (because we smooth for velocity filtering)

%{
% cut each elapsed time in half and assign it to each coordinate value.
% that way each coordinate value gets half of the time before it and half
% of the time after it
betweens=timestamps-(dt/2); % not sure i want to do this...
betweens(end+1)=coorddata(end,1);
elapsed=diff(betweens);
%}
elapsed=dt;
elapsed(elapsed>timejumpdelay)=nan;

if isempty(VelData)
% this is the smoothing algorithm to get the velocity data
windowspan=round(velsmooth/median(dt))*5; % how many seconds to smooth the data across?
coorddata2=[SmoothMat2(coorddata(:,2),[1 windowspan*4], windowspan),...
    SmoothMat2(coorddata(:,3),[1 windowspan*4], windowspan)];


% this is in cm
displacement=sqrt(diff(coorddata2(:,1)).^2+diff(coorddata2(:,2)).^2)/PixPerCM;
%this is a non smoothed displacement, we may want to bin this to remove the
%small jitters that occur
% any displacement where the time is too far, we have to nan out

% nan out long lags so they dont contribute to the smoothed average
displacement(elapsed>timejumpdelay)=nan; 
% for smoothing the velocity data, we grab all the coords that are within
% our timewindow average them and thats our smoothed velocity.
velocity=nan(length(displacement),1);

% this takes fucking forever, so lets just do a 5 bin arma
parfor i=1:length(velocity)
    velocity(i)=nansum(displacement(abs(timestamps-timestamps(i))<=velsmooth))/...
        nansum(elapsed(abs(timestamps-timestamps(i))<=velsmooth));
end

%velocity=smooth(displacement./dt,vfactor,'moving');
% there will be some high velocity timestamps because the strobe jumps a
% bit, generally about 60 jumps between 2 and 10 seconds
else
    velocity=VelData;
end



velcoords = (cat(2,timestamps,xycoords,velocity,elapsed));


%coords now with velocities and latencies
% cols as follows
% ts  |  xx   |  yy  |  vel  |   dt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% now threshold based on velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
velcoords = velcoords(velcoords(:,4) >= minvelocity,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% disregard time jumps%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% no we we have a ton of snippits of fast moving data.  We want to make
% sure that we dont take spikes that are outside of those data, and we want
% to make sure we dont count elapsed time outside of those data
% we also dont want to use high velocity epochs that are super quick,
% because thoes may be just head pans
% this is in seconds
SpikeTimes=unit.ts(:,1);

timejumps = [];
timejumps = cat(1,timejumps,find(diff(velcoords(:,1)) > timejumpdelay));
% data isnt interpolated, so we will have to use a timejump window of
% around 10, not .334 (which is strobe rate)
binspikes = [];
session_duration = 0;
% log session duration and pull spikes that are ONLY within our epochs,
% epochs, that is, that have ocntinuous tracking (diff is less than 10
% second
% this only remove spikes, this doesnt modify any coordinates
if length(timejumps) > 0
    for i = 1:length(timejumps)+1
        if i == 1 % if its the first epoch
            % remove ts that will be assigned earlier than ts1 and after
            % the ts before your pause
            binspikes = cat(1,binspikes,SpikeTimes(SpikeTimes >= velcoords(1) & SpikeTimes <=velcoords(timejumps(i),1)));
            session_duration = timestamps(timejumps(i))-timestamps(1);
        elseif  i == length(timejumps)+1 % if its the last epoch
            % remove the ts that will be assigned to the first ts of the
            % last epoch
            binspikes = cat(1,binspikes,SpikeTimes(SpikeTimes >= velcoords((timejumps(i-1)+2),1) & SpikeTimes <= velcoords(end,1)));
            session_duration = session_duration+(timestamps(length(timestamps))-(timestamps((timejumps(i-1)+2))));
        else % middle epochs
            % reomove spikes that will be assigned to the ts that border
            % the pause (each will be assigned a large lapsed time
            binspikes = cat(1,binspikes,SpikeTimes(SpikeTimes >= velcoords((timejumps(i-1)+2),1) & SpikeTimes <= velcoords(timejumps(i),1)));
            session_duration = session_duration+velcoords(timejumps(i),1)-velcoords((timejumps(i-1)+2),1);
        end
    end
    newSpikeTimes = binspikes;
    
else
    newSpikeTimes = SpikeTimes(SpikeTimes >= velcoords(1) & SpikeTimes <= velcoords(end,1));
    session_duration =  velcoords(end,1)-velcoords(1);
end

% after taking out dead space, what do we have left?
tossedspikes=length(SpikeTimes)-length(newSpikeTimes);
if suppress<2
    fprintf('tossing %d of a total %d spikes \n', tossedspikes, length(unit.ts));
end
% just an error check
if find(diff(sort(SpikeTimes))) == 0 %double check to make sure timestamps aren't duplicated
    disp 'Error: Repeated Spike Time'
    
end


if length(newSpikeTimes)<spikethreshold
    % show warnings only if you want
    if suppress<2
        fprintf('cell has fewer than %d spikes \n', spikethreshold);
    end
    ratemap=[];
    %g=figure;
else
    average_fr=length(newSpikeTimes)/session_duration;


    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Now plot the spikes and trajectory ts%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % now we've ditched any spikes outside of 0.3 seconds from the closest
    % video timestamp
    % cellspikes are now ts, x, y, vel, elapsed time
    cellspikes = [newSpikeTimes interp1(velcoords(:,1),velcoords(:,2:end),newSpikeTimes,'nearest')];
    
    topplot.spk=cellspikes(:,2:3);
    topplot.occup=velcoords(:,2:3);
    %%%%%%%%%%%%%%%% save top plot data %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is before binning
    
    
    
    % top plot will be the spikes and occupancy
    if suppress==0
        g=figure;
        if onlyheat==0
        m=subplot(2,1,1);
        
        plot(velcoords(:,2), velcoords(:,3), 'k*');
        hold on;
        %rat path for entire session, regardless of velocity
        
        plot(cellspikes(:,2), cellspikes(:,3), 'r.');
        hold off
        %axis image
        axis off
        set(gca, 'Box', 'on')
        end
    end
    
    
    topplot.occ(:,1)=velcoords(:,2); topplot.occ(:,2)=velcoords(:,3);
    
    %title('Entire Session','fontsize',14,'FontName', 'Arial');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get resolution for bottom graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % realcoords is if you have a 480 to 640 frame size (usually its
    % higher)
    if realcoords==1
         XX=0:(PixPerCM*CMPerBin):1600; % probably need to make this bigger
         YY=0:(PixPerCM*CMPerBin):1600;
    else
        % this is if you dont knoiw your whole window size
        XX=min(velcoords(:,2)):PixPerCM*CMPerBin:max(velcoords(:,2));
        YY=min(velcoords(:,3)):(PixPerCM*CMPerBin):max(velcoords(:,3));
    end
    % now resize our coordinates, fit each coordinate into a bin set by our
    % factor
    for i=1:length(velcoords)
        [~,xindices(i)]=min(abs(XX-velcoords(i,2)));
        [~,yindices(i)]=min(abs(YY-velcoords(i,3)));
    end
    
    % Sometimes nan coordinates carry throguh and you need to delete these
    nanoverx=isnan(velcoords(:,2)); nanovery=isnan(velcoords(:,3));
    velcoords(:,2)=xindices; velcoords(:,3)=yindices;
    velcoords(nanoverx,2)=nan; velcoords(nanovery,3)=nan; 
    % now velcoords are in place map bins
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%% Accumulate bin counts occup/spk %%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % the accumarray option
    % accumarray across all bins because you want to account for zero visit
    % bins, this accumulates 5th column which is elapsed time (dt)
    noccup=accumarray([xindices',yindices'],velcoords(:,5),[length(XX),length(YY)]);

    % cast the spikes into the same timestamps but update the coordinates

    cellspikes = [newSpikeTimes interp1(velcoords(:,1),velcoords(:,2:end),newSpikeTimes,'nearest')];
    % hist3 is the same as accumarray
    newedges={1:length(XX),1:length(YY)};
    [nspikes,~]=hist3(cellspikes(:,2:3),'Edges',newedges);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Now filter pixels based on visit counts and time spent %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % calculating visits:
    I_X=zeros(1,length(velcoords(:,2)));
    I_Y=I_X;
    for j=1:length(velcoords(:,2))
        [~,I_X(j)]=min(abs(velcoords(j,2)-newedges{1}));
        [~,I_Y(j)]=min(abs(velcoords(j,3)-newedges{2}));
    end
    dx=diff([I_X(1)-1,I_X]); %this is to add a value for the first datapoint
    dy=diff([I_Y(1)-1,I_Y]);
    % i.e. if the coordinates are the same, its the same visit..
    % generally you should pad this...
    % not sure why he doesnt just take abs of each of these...
    newvisitstimes=-(dy == 0).*(dx==0)+1;
    
    
    [nvisits,~]=hist3([velcoords(newvisitstimes>0,2),velcoords(newvisitstimes>0,3)],'Edges',newedges);

    
    % nix pixels that werent occupied long enough
    goodpix=noccup>=mintimespent;
    % nix pixels that werent visited enough times
    Valid=(nvisits>minvisits).*goodpix;
    % and fills the holes in that box
    Valid=imfill(Valid);

    % valid is now what we will use to plot occupancy, and to remove pixels
    % after the grid is smoothed
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% now i smooth %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %{
    There are two options here, first to convolve separately, which
     is how the mosers did it back in their 05? paper
    the other way is to convolve after the fact. which may be more
    appropriate because its most likely that outlier pixels for spikerate
    will also be outlier pixels for occupancy. I have not found that an
    outlier in occupancy has an outlier in rate, but if thats hte case then
    convolve separately

    
    % we will have to smooth occupancy bur remove bad pixels;

    
    %%%%%% heres the code to smooth before you divide spikes by occupancy
    
    %smoothoccup = SmoothMat(noccup, [Gausswin,Gausswin], Gaussdev);

    %smoothspikes = SmoothMat(nspikes, [Gausswin,Gausswin], Gaussdev);
    
    
    % now get ratemap
    %ratemap=smoothspikes./smoothoccup;
    
    %}    
    
    
    
    % take only valid pixels
    nspikes=nspikes.*Valid; noccup=noccup.*Valid;
    ratemap=nspikes./noccup;
    
    % Save out raw ratemap
    rawratemap=ratemap.*Valid; rawratemap(isnan(rawratemap))=0;
    
    % smooth the map here
    ratemap=SmoothMat2(ratemap,[Gausswin,Gausswin], Gaussdev);
    
    % and fill in where there is no occupancy (triple check
    ratemap=ratemap.*Valid;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Clean up the image %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %i have to transpose all matrices now so their dimensions agree
    alphadata = (Valid' == 0);
    ratemap=ratemap';
    % this is an output so it should agree too
    noccup2=noccup';
    
    % create a white background
    white = cat(3, ones( size(ratemap)), ones( size(ratemap)), ones( size(ratemap)));
    gray = white.*.8;
    floorrate=.01;
    nospikes=(ratemap<=floorrate);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %%%%%  FIRING RATE LIMITS
    %%%%%%%%%%%%%%%%%%%%%
    
    % find all pixel rates
    rates=reshape(ratemap,1,numel(ratemap));
    
    % parse down to only real rates
    realrates=rates(~isnan(rates) & ~isinf(rates));
    
    % now pad the top 97% of pixels
    if RealLim
        if isempty(clim)
            clim=max(max(realrates));
        end
    else
        if isempty(clim)
            [f,x]=ecdf(realrates);
            clim=median(x(f>.97));
        end
    end
    
    % now any high firing rate is turned to the clim (artificial max)
    for w=1:numel(ratemap)
        if ratemap(w)>=clim
            ratemap(w)=clim;
        end
    end
    
    %
    
    colormap_fr =linspace(0,clim,255);
    [~,colormap]=histc(ratemap,colormap_fr);
    rgbcolormap=ind2rgb(colormap,jet(256));
    
    if grayout
        finalcolormap1=immerge(rgbcolormap,gray,nospikes);
        finalcolormap=immerge(finalcolormap1,white,alphadata);
    else
        finalcolormap=immerge(rgbcolormap,white,alphadata);
    end
    if suppress==0
        if onlyheat==0
            n=subplot(2,1,2,'align');
        end
        q=image(finalcolormap);
        set(gcf,'Colormap',jet(256));
        axis off
        set(gca,'YDir','normal')
        FRticks=linspace(0,1,6);
        FRinterval=round(linspace(0,clim,6),1);
        f=colorbar; set(f,'YTick',FRticks,'YTickLabel',FRinterval);
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % now cosmetic stuff%%%
        %%%%%%%%%%%%%%%%%%%%%%%
        if onlyheat==0
        % to get the plots to line up;
        bottompos=get(n,'Position');
        figwidth=bottompos(3);
        toppos=get(m,'Position');
        toppos(3)=figwidth;
        set(m,'Position',toppos);
        % matlab has a broken alphadata function, so we have to immerge our
        % colormaps.  Maybe in future versions of matlab this wont be a problem
        %set(q,'AlphaData',noccup');
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% add bwboundary to create contiguous place fields%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        %axis off
        %axis image
        
        
        % this puts white in the background
        try
            title({[' unit ' num2str(unit.units) ' ' EpochName],['mean firing rate = ' num2str(average_fr) 'Hz']} );
        end
    end
end
end



% just a function to plot shit
function plotitout



smoothx=movsum(coorddata(:,2),12,'Endpoints','fill')/12;
smoothy=movsum(coorddata(:,3),12,'Endpoints','fill')/12;
displacement2=sqrt(diff(smoothx).^2+diff(smoothy).^2)/PixPerCM;
velocity2=smooth(displacement2./dt,vfactor,'moving');

sp=subplot(3,1,1); plot(timestamps(1:10000),xycoords(1:10000,1));
hold on;
plot(timestamps(1:10000),xycoords(1:10000,2))
sp(2)=subplot(3,1,2);
plot(timestamps(1:10000),displacement(1:10000));
hold on
plot(timestamps(1:10000),displacement2(1:10000));
sp(3)=subplot(3,1,3);
plot(timestamps(1:10000),velocity(1:10000));
hold on;
plot(timestamps(1:10000),velocity2(1:10000));
linkaxes(sp,'x');

end

