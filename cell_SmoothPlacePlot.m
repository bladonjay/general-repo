function [ratemap,topplot,finalcolormap,noccup2,cellspikes,velcoords,nspikes,noccup] = cell_SmoothPlacePlot(session,unit,varargin)
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
% this function requires the following functions:
%

%
% JHB 9-5-14

p=inputParser;
addOptional(p,'Gausswin',40); % huge cause it drops off real quick
addOptional(p,'Gaussdev',1);
addOptional(p,'Factor',6); % meaning pixels per bin
addOptional(p,'minvelocity',20); % in pixels per second
addOptional(p,'mintimespent',.05);  % or about 5 timestamps
addOptional(p,'minvisits',2);
addOptional(p,'velsmooth',.5); % in seconds
addOptional(p,'spikethreshold',30);
addOptional(p,'clim',[]);
addOptional(p,'EpochName','Whole Session');
addOptional(p,'RealLim',0);
addOptional(p,'Grayout',0);
addOptional(p,'Timejumpdelay',1);
addOptional(p,'suppress',0);
addOptional(p,'onlyheat',0);
addOptional(p,'realcoords',0);
parse(p,varargin{:});


%%%%% these are plotting parameters %%%%%%%
% fewest spikes the cell has to have to plot
spikethreshold=p.Results.spikethreshold;
% size of the gaussian kernel
Gausswin=p.Results.Gausswin;
% STD of gaussian kernel
Gaussdev=p.Results.Gaussdev;
% this is the pixel size (play with it)
Factor=p.Results.Factor*6.667;
% if you want to scale the rateplot based on a set firing rate max, e.g.
% when you are producing paired plots
clim=p.Results.clim;
%use max firing rate or max reasonable?
RealLim=p.Results.RealLim;
% gray out no spikes?
grayout=p.Results.Grayout;
% Name the epoch for the figure
EpochName=p.Results.EpochName;
% suppress plot
suppress=p.Results.suppress;
% only heat plot
onlyheat=p.Results.onlyheat;

%%%%% these are cutoff parameters %%%%%%%%
% velocity cutoff to remove any stationary time
minvelocity=p.Results.minvelocity;
% min time spent in each pixel
mintimespent=p.Results.mintimespent;
% min visits to the pixel
minvisits=p.Results.minvisits;
% time in seconds to smooth the velocity
velsmooth=p.Results.velsmooth;
% time delay in seconds that is too long to interpolate
timejumpdelay=p.Results.Timejumpdelay;

% real coordinates so that you can compare pixels
realcoords=p.Results.realcoords;



% output variables;
ratemap=[];
rawratemap=[];

finalcolormap=[];
noccup2=[];
cellspikes=[];
velcoords=[];
nspikes=[];
noccup=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get velocity and make velocity filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% remove the first coordinate
tcoord = session.edit_coords((2:end),1);
tempcoords = session.edit_coords((2:end),2:3);


%difference between consecutive coords, but this is for before each
%timestamp. e.g. dt(1) is the difference between tcoord(0) and tcoord(1),
%or the elapsed time BEFORE each timestamp
dt = diff(session.edit_coords(:,1));
% this is the time window that continuous time data will be assigned that
% tcoord

% this is if time were continuous, the elapsed is the length of time
% assigned to that ts


displacement=sqrt(diff(session.edit_coords(:,2)).^2+diff(session.edit_coords(:,3)).^2)/median(dt);
% remove if youre only moving one or two pixels over
%this is a non smoothed displacement, we may want to bin this to remove the
%small jitters that occur, lets use a moving average

% vfactor is the amount of time to smooth
% probably should discretize this dataset so we can bin low velocities to
% zero and get a decent spread of velocities * function discretize*
% smooth factor for velocity must be odd
vfactor=round((velsmooth/median(dt))/2)*2-1;
rawvelocity=smooth(displacement,vfactor,'moving'); % mnaybe a gaussian window
velocity=(rawvelocity.*dt)* 30; % in pixels per second
% there will be some high velocity timestamps because the strobe jumps a
% bit, generally about 60 jumps between 2 and 10 se conds
velcoords = (cat(2,tcoord,tempcoords(:,1),tempcoords(:,2),velocity,dt));
%coords now with velocities and latencies
% cols as follows
% ts  |  xx   |  yy  |  vel  |   dt

% we can kalman filter velocity if we want
SpikeTimes=unit.ts(:,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% disregard time jumps %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you want this done right you will have to look at a histogram of time
% jumps, this is easily done by plotting or histogramming diff(timestamps)

% this is in seconds

% timejumps is the index of the coordinate BEFORE the delay
timejumps = [];
timejumps = cat(1,timejumps,find(diff(tcoord) > timejumpdelay));
% data isnt interpolated, so we will have to use a timejump window of
% around 10, not .334 (which is strobe rate)
binspikes = [];
session_duration = 0;
% log session duration and pull spikes that are ONLY within our epochs,
% epochs, that is, that have ocntinuous tracking (diff is less than 10
% second
% this only remove spikes, this doesnt modify any coordinates
if length(timejumps) >= 1
    for i = 1:length(timejumps)+1
        if i == 1 % if its the first epoch
            % remove ts that will be assigned earlier than ts1 and after
            % the ts before your pause
            binspikes = cat(1,binspikes,SpikeTimes(SpikeTimes >= tcoord(1) & SpikeTimes <= tcoord(timejumps(i))));
            session_duration = tcoord(timejumps(i))-tcoord(1);
        elseif  i == length(timejumps)+1 % if its the last epoch
            % remove the ts that will be assigned to the first ts of the
            % last epoch
            binspikes = cat(1,binspikes,SpikeTimes(SpikeTimes >= tcoord((timejumps(i-1)+2)) & SpikeTimes <= tcoord(end)));
            session_duration = session_duration+(tcoord(length(tcoord))-(tcoord((timejumps(i-1)+2))));
        else % middle epochs
            % reomove spikes that will be assigned to the ts that border
            % the pause (each will be assigned a large lapsed time
            binspikes = cat(1,binspikes,SpikeTimes(SpikeTimes >= tcoord((timejumps(i-1)+2)) & SpikeTimes <= tcoord(timejumps(i))));
            session_duration = session_duration+(tcoord(timejumps(i))-(tcoord((timejumps(i-1)+2))));
        end
    end
    newSpikeTimes = binspikes;
    
else
    newSpikeTimes = SpikeTimes(SpikeTimes >= tcoord(1) & SpikeTimes <= tcoord(end));
    session_duration =  tcoord(length(tcoord)) - tcoord(1);
end

%%%%%%%% reconstruct our valid spikes %%%%%%%%%%%

% our new cell spike times
cellspikes=newSpikeTimes;

% after taking out dead space and scrubbing, what do we have left?
tossedspikes=length(SpikeTimes)-length(cellspikes);
if suppress<2
    fprintf('tossing %d of a total %d spikes \n', tossedspikes, length(unit.ts));
end
% just an error check
if find(diff(SpikeTimes)) == 0 %double check to make sure timestamps aren't duplicated
    disp 'Error: Repeated Spike Time'
end




% run the spikes, but only if the cell has enough spikes
if length(cellspikes)<spikethreshold
    % show warnings only if you want
    if suppress<2
        %fprintf('cell %s has fewer than %d spikes \n', unit.name, spikethreshold);
        fprintf('cell has fewer than %d spikes \n', spikethreshold);
    end
    topplot.occ=[0 0];
    topplot.spk=[0 0];
    
    ratemap=[];
    %g=figure;
    return
else
    
    % if we have enough spikes lets roll
    
    %%%%%%%%%%%% Find pixels for each spike and scrub for long delays
    
    % cell spikes is as follows:
    % ts | x coord | y coord | velocity | dt
    cellspikes = interp1(velcoords(:,1),velcoords(:,1:end),newSpikeTimes,'nearest');
    cellspikes = cellspikes(cellspikes(:,4) > minvelocity,:);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Now plot the spikes and occupancy ts %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % take out when the velocity is actually zero from both occupancy and
    % spikes
    velcoords=velcoords(velcoords(:,4)>minvelocity,:);
    % top plot will be the spikes and occupancy
    if suppress==0
        g=figure;
        if onlyheat==0
            m=subplot(2,1,1);
            
            plot(velcoords(:,2), velcoords(:,3), 'k-');
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
    topplot.spk(:,1)=cellspikes(:,2); topplot.spk(:,2)=cellspikes(:,3);
    
    
    %title('Entire Session','fontsize',14,'FontName', 'Arial');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get resolution for bottom graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%% REDO THIS, USE HISTC as below to get your pixels right****
    % decimate our spatial sampling based on pix/cm
    rangex=linspace(min(velcoords(:,2)),max(velcoords(:,2)),1200/(Factor));
    % now histc it
    [~,indx]=histc(velcoords(:,2),rangex);
    % fill in real values
    binnedx=rangex(indx);
    % indices will be valid, youre just interval scaling it now
    xpixels=indx-min(indx)+1;
    
    % same for y
    rangey=linspace(min(velcoords(:,3)),max(velcoords(:,3)),1200/(Factor));
    [~,indy]=histc(velcoords(:,3),rangey);
    binnedy=rangey(indy);
    ypixels=indy-min(indy)+1;
    %     [bins,bininds]=histc(topplot.occ(:,1),binrange);
    %     occ2(:,1)=binrange(bininds);
    
    
    
    %
    %     % factor is literally centimeters per pixel
    %     % you'll have to check this for the camera, but in our rooms, a pixel
    %     % is about 6 cm, so
    %     if realcoords==1
    %         XX=0:(Factor):1200;
    %         YY=0:(Factor):1200;
    %     else
    %         % here is if youd like to start at 0 for your first coordinate,
    %         % e.g. where you decimate the dataset
    %
    %         XX=min(velcoords(:,2)):(Factor):max(velcoords(:,2));
    %         YY=min(velcoords(:,3)):(Factor):max(velcoords(:,3));
    %     end
    %     % now resize our coordinates, fit each coordinate into a bin set by our
    %     % factor
    %     for i=1:length(velcoords)
    %         [~,xindices(i)]=min(abs(XX-velcoords(i,2)));
    %         [~,yindices(i)]=min(abs(YY-velcoords(i,3)));
    %     end
    
    
    
    
    
    
    
    
    
    % now snap coords to pixels to get occupancy
    
    noccup=accumarray([xpixels ypixels],velcoords(:,5),[max(xpixels),max(ypixels)]);
    % this sums up all the dt times for each coordinate
    % get noccup to be real time;
    
    
    
    % now snap spikes to their appropriate pixel
    
    velcoords(:,2)=xpixels; velcoords(:,3)=ypixels;
    
    newedges={1:max(xpixels),1:max(ypixels)};
    
    
    % now rebuild our cellspikes for Coords
    cellspikesC = interp1(velcoords(:,1),velcoords(:,1:end),newSpikeTimes,'nearest');
    % same filter
    cellspikesC = cellspikesC(cellspikes(:,4) > minvelocity,:);
    
    % getting the number of spikes in each bin
    [nspikes,~]=hist3(cellspikesC(:,2:3),'Edges',newedges);
    
    
    %%%%%% noccup is a grid of occupancy (in seconds), and nspikes is a grid of spike
    %%%%%% locations %%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Now filter based on pixels %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
    newvisitstimes=double(abs(dy) > 0 | abs(dx)>0);
    
    
    [nvisits,~]=hist3([velcoords(newvisitstimes>0,2),velcoords(newvisitstimes>0,3)],'Edges',newedges);
    
    
    % now occupation
    
    
    % now pixels that have really low visits and really high firing rates;
    highfr=reshape(nspikes,1,numel(nspikes));
    highfr=highfr(highfr~=0); % to only use pixels with spikes (theres alot of dead space)
    % these are the pixels where the firing rate was top 5%
    %[f,x]=ecdf(highfr);
    %highrate=mean(x(f>.95)); % find the average of all the high pixels
    
    % nix pixels that werent occupied
    badpix=(noccup<mintimespent);
    goodpix=~badpix;
    
    
    
    
    % nix pixels that werent visited enough times
    Valid=(nvisits>minvisits).*goodpix;
    % and fills the holes in that box
    Valid=imfill(Valid);
    
    
    
    % valid is now what we will use to plot occupancy, and to remove pixels
    % after the grid is smoothed
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% now i smooth %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % There are two options here, first to convolve separately, which
    %  is how the mosers did it back in their 05? paper
    % the other way is to convolve after the fact. which may be more
    % appropriate because its most likely that outlier pixels for spikerate
    % will also be outlier pixels for occupancy. I have not found that an
    % outlier in occupancy has an outlier in rate, but if thats hte case then
    % convolve separately
    
    
    % we will have to smooth occupancy bur remove bad pixels;
    nspikes=nspikes.*Valid;
    noccup=noccup.*Valid;
    
    %%%%%% heres the code to smooth before you divide spikes by occupancy
    
    %smoothoccup = SmoothMat(noccup, [Gausswin,Gausswin], Gaussdev);
    
    % we can smooth all spikes because they all will be valid
    %smoothspikes = SmoothMat(nspikes, [Gausswin,Gausswin], Gaussdev);
    
    
    % now get ratemap
    %ratemap=smoothspikes./smoothoccup;
    
    ratemap=nspikes./noccup;
    rawratemap=ratemap.*Valid;
    rawratemap(isnan(rawratemap))=0;
    %ratemap(isnan(ratemap))=0; % in case youre not using nanconv, this will
    %smooth down all edges
    ratemap=SmoothMat2(ratemap,[Gausswin,Gausswin], Gaussdev);
    
    % and fill in where there is no occupancy
    ratemap=ratemap.*Valid;
    
    
    
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
        if isempty(clim);
            clim=max(max(realrates));
        end
    else
        if isempty(clim);
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
        q=image(finalcolormap); %,'AlphaData',noccup');
        set(gcf,'Colormap',jet(256));
        axis off
        set(gca,'YDir','normal')
        %colorbar
        FRinterval=linspace(0,clim,6);
        %f=colorbar; set(f,'YTickLabel',FRinterval(2:end));
        title(sprintf('max %.2f',clim));
        
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
        

    end
end
end

