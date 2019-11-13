function [ratemap,topplot,finalcolormap,noccup2,cellspikes,velcoords,rawdata] = SmoothPlacePlot2_OLD(tracking,unit,varargin)
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

%
% JHB 9-5-14

p=inputParser;
addOptional(p,'Gausswin',20);
addOptional(p,'Gaussdev',1);
addOptional(p,'Factor',12);
addOptional(p,'minvelocity',0.1);
addOptional(p,'mintimespent',.2);  % or about 5 timestamps
addOptional(p,'minvisits',1);
addOptional(p,'velsmooth',.2); % in seconds
addOptional(p,'spikethreshold',50);
addOptional(p,'clim',[]);
addOptional(p,'EpochName','Whole Session');
addOptional(p,'RealLim',1);
addOptional(p,'Grayout',0);
addOptional(p,'Timejumpdelay',1);
addOptional(p,'suppress',0);
addOptional(p,'onlyheat',0);
addOptional(p,'realcoords',0);
addOptional(p,'smoothfirst',0);
parse(p,varargin{:});


%%%%% these are plotting parameters %%%%%%%
% fewest spikes the cell has to have to plot
spikethreshold=p.Results.spikethreshold;
% size of the gaussian kernel
Gausswin=p.Results.Gausswin;
% STD of gaussian kernel
Gaussdev=p.Results.Gaussdev;
% this is the pixel size (play with it)
Factor=p.Results.Factor;
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
% smooth before ratemap
smoothfirst = p.Results.smoothfirst;

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
topplot=[];
finalcolormap=[];
noccup2=[];
cellspikes=[];
velcoords=[];
rawdata.ratemap=nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get velocity and make velocity filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
strobe=tracking(:,1);
dt = [diff(strobe); .033]; %difference between consecutive coords

velcoords = tracking(:,2:3);
% this is if time were continuous, the elapsed is the length of time
% assigned to that ts

displacement=[sqrt(diff(velcoords(:,1)).^2+diff(velcoords(:,2)).^2)*Factor; 0];
%this is a non smoothed displacement, we may want to bin this to remove the
%small jitters that occur


% smooth factor for velocity must be odd
vfactor=round((velsmooth/min(dt))/2)*2-1;
velocity=smooth(displacement.*dt,vfactor,'moving');
% there will be some high velocity timestamps because the strobe jumps a
% bit, generally about 60 jumps between 2 and 10 seconds
%coords now with velocities and latencies
% cols as follows
% ts  |  xx   |  yy  |  vel  |   dt
velcoords=[strobe velcoords velocity dt];
% we can kalman filter velocity if we want
SpikeTimes=unit.ts;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% disregard time jumps%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you want this done right you will have to look at a histogram of time
% jumps, this is easily done by plotting or histogramming diff(timestamps)

% this is in seconds


timejumps = [];
% the indices of the time jumps
timejumps = [0; cat(1,timejumps,find(dt > timejumpdelay)); length(dt)];
% data isnt interpolated, so we will have to use a timejump window of
% around 10, not .0334 (which is strobe rate)
binspikes = [];
session_duration = 0;
% log session duration and pull spikes that are ONLY within our epochs,
% epochs, that is, that have ocntinuous tracking (diff is less than 10
% second
% this only remove spikes, this doesnt mod
for i = 1:length(timejumps)-1
    
    binspikes = cat(1,binspikes,SpikeTimes(SpikeTimes >= strobe(timejumps(i)+1) & SpikeTimes <= strobe(timejumps(i+1))));
    session_duration = session_duration+(strobe(timejumps(i+1))-(strobe((timejumps(i)+1))));
end

newSpikeTimes = binspikes;


% after taking out dead space, what do we have left?
tossedspikes=length(SpikeTimes)-length(newSpikeTimes);
if suppress<2
    fprintf('tossing %d of a total %d spikes \n', tossedspikes, length(unit.ts));
end
% just an error check
if find(diff(SpikeTimes)) == 0 %double check to make sure timestamps aren't duplicated
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
    average_fr=length(SpikeTimes)/session_duration;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% now take out timestamps where the dwell time is too long
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    % find out the new ts, X,Y, velocity, and duration for each spike, because the
    % spiketimes dont exactly match up with the time coordinates, you can
    % interpolate, but this just takes the closest values
    
    

    % cell spikes is as follows:
    % ts | x coord | y coord | velocity | dt
    
    % looks like this is the best place to pull out the velocity to see if
    % there is any velocity coding
    
    % if the cell has too few spikes, break the functio
    
    %%%%%%%% here you will end with two four column matrices, one for the cell
    %%%%%%%% and one for the rat.
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Now plot the spikes and occupancy ts%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % this is before the resize, cellspikes will be overwritten
    cellspikes = interp1(velcoords(:,1),velcoords(:,2:end),newSpikeTimes,'nearest');
    cellspikes=[newSpikeTimes cellspikes]; % this just tacks on the timestamps to the rest of the data
    % now filter based on velocity
    cellspikes = cellspikes(cellspikes(:,4) >minvelocity,:);
    
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
    topplot.spk(:,1)=cellspikes(:,2); topplot.spk(:,2)=cellspikes(:,3);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get resolution for bottom graph
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % factor is literally centimeters per pixel
    % basically you can divide your numbers by a factor and round
    if realcoords==1
         XX=0:(Factor):1200;
         YY=0:(Factor):1200;
    else
        % here is if youd like to start at 0 for your first coordinate
        XX=min(velcoords(:,2)):(Factor):max(velcoords(:,2));
        YY=min(velcoords(:,3)):(Factor):max(velcoords(:,3));
        velcoords(:,2)=velcoords(:,2)-min(velcoords(:,2));
        velcoords(:,3)=velcoords(:,3)-min(velcoords(:,3));
    end
    % now resize our coordinates, fit each coordinate into a bin set by our
    % factor
    
    [~,xindices]=histc(velcoords(:,2),XX);
    [~,yindices]=histc(velcoords(:,3),YY);
    
    velcoords(:,2)=xindices; velcoords(:,3)=yindices;
    
    % now interpolate the xy coords for our spikes
    cellspikes = interp1(velcoords(:,1),velcoords(:,2:end),newSpikeTimes,'nearest');
    % tack on their ts
    cellspikes=[newSpikeTimes cellspikes];
    
    
    % now filter based on minimum speed, same as above, the crucial bit
    % here is that we dont recalculate the dwell time, cause we're
    % scrubbing data
    cellspikes = cellspikes(cellspikes(:,4) >minvelocity,:);
    % now turn it into a grid
    
    % getting the occupancy in number of times the animal visited the bin
    noccup=accumarray([xindices yindices],velcoords(:,5),[length(XX),length(YY)]);
    % this sums up all the dt times for each coordinate to get
    % noccup to be real time
    
    newedges={1:length(XX),1:length(YY)};
    
    % getting the number of spikes in each bin
    [nspikes,~]=hist3(cellspikes(:,2:3),'Edges',newedges);
    
    %%%%%%%%%%%%%%%% HERE WE HAVE raw occupancy and raw spike counts per
    %%%%%%%%%%%%%%%% bin

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Now filter based on pixels %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % calculating visits:
    % preallocate space
    
    dx=diff([xindices(1)-1; xindices]); %this is to add a value for the first datapoint
    dy=diff([yindices(1)-1; yindices]);
    % i.e. if the coordinates are the same, its the same visit..
    % generally you should pad this...
    % not sure why he doesnt just take abs of each of these...
    newvisitstimes=-(dy == 0).*(dx==0)+1;
    [nvisits,~]=hist3([velcoords(newvisitstimes>0,2),velcoords(newvisitstimes>0,3)],'Edges',newedges);
    
    
    % now occupation
    % now pixels that have really low visits and really high firing rates;
    % first find the pixels with really high firing rates
    %highpix=nspikes-nanmean(nspikes(:))>2*stdev(nspikes(:))
    
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
    % is how the mosers did it back in their 05 paper
    % the other way is to convolve after the fact. Its most likely that outlier pixels for spikerate
    % will also be outlier pixels for occupancy. I have not found that an
    % outlier in occupancy isnt also an outlier in rate, but if thats the case then
    % convolve separately
    
    
    % we will have to smooth occupancy bur remove bad pixels;
    nspikes=nspikes.*Valid;
    noccup=noccup.*Valid;
    
    rawdata.nspikes=nspikes; rawdata.noccup=noccup;
    %%%%%% heres the code to smooth before you divide spikes by occupancy
    if smoothfirst
        
        smoothoccup = SmoothMat(noccup, [Gausswin,Gausswin], Gaussdev);
        
        % we can smooth all spikes because they all will be valid
        smoothspikes = SmoothMat(nspikes, [Gausswin,Gausswin], Gaussdev);
        
        
        % now get ratemap
        ratemap=smoothspikes./smoothoccup;
    else
        ratemap=nspikes./noccup;
        rawratemap=ratemap.*Valid;
        rawratemap(isnan(rawratemap))=0;
        rawdata.ratemap=rawratemap;
        %ratemap(isnan(ratemap))=0; % in case youre not using nanconv, this will
        %smooth down all edges
        ratemap=SmoothMat2(ratemap,[Gausswin,Gausswin], Gaussdev);
    end
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
    % if we want the real max pixel and we havent designated a max, use it
    if RealLim
        if isempty(clim)
            clim=max(max(realrates));
        end
        % if we want to pad, take anything bigger than our clim and pad it
        % to max (if no clim, find max pixel)
    else
        if isempty(clim)
            [f,x]=ecdf(realrates(realrates>0));
            clim=mean(x(f>.97));
        end
    end
    
    % now any high firing rate is turned to the clim (artificial max)
    ratemap(ratemap>clim)=clim;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Get Skaggs Information %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try
    goodpix=noccup2>0;
    pixprob=linearize(noccup2(goodpix))./sum(sum(noccup2(goodpix)));
    SpatialInfo=Skaggs_basic(pixprob,linearize(ratemap(goodpix)),...
        nanmean(nanmean(ratemap(goodpix))));
    catch
        SpatialInfo=nan;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Get place field size %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PlaceMax=clim*.8; hotnes=[];
    try
        patchmat=ratemap>PlaceMax;
        details = regionprops(patchmat,'PixelIdxList','Area');
        for ds=1:length(details)
            hotness(ds)=nanmean(ratemap(details(ds).PixelIdxList));
        end
        [~,bigfield]=max(hotness);
        PlaceFieldSize=details(bigfield).Area;
        % no get hottest region and get its pixel size
    catch
        PlaceFieldSize=nan;
    end
    rawdata.PlaceFieldSize=PlaceFieldSize;
    rawdata.SpatialInfo=SpatialInfo;
    rawdata.MaxRate=clim;
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
        image(finalcolormap);
        text(.2,-.1,sprintf('Max FR: %.2f, PF Size : %.2f, Bits/Spike: %.2f',clim, PlaceFieldSize,SpatialInfo),...
            'units','normalized');
        set(gcf,'Colormap',jet(256));
        axis off
        set(gca,'YDir','normal');
        f=colorbar; set(f,'YTickLabel',round(colormap_fr(1:25:end),2));
        
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
        
    end
end
end

