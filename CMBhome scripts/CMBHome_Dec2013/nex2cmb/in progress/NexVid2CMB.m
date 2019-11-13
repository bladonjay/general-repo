function [viddata, nex] = NexVid2CMB(nexfile,varargin)
%Imports any and all information from a cineplex NEX file into the
%CMBHome data class.
%NexVid2CMB = (nex file, varargin)
% varargin can be a cell struct (containing a list of events to keep), a
% logical (1 for batch mode), and 'nvt' or 'nlx' if its a neuralynx file

% JHB 6-10-14

eventstokeep = {};
batch = false;
nlx = false;


for a = 1:numel(varargin)
    if(iscellstr(varargin{a}))
        % 'eventstokeep' dictates which events from the NEX file are copied into
        % the CMBHOME object. If empty, all events are kept. To exclude all events,
        % include one element in 'eventstokeep' that won't match any other event,
        % for example: eventstokeep = {''};
        eventstokeep = lower(varargin{a});
    elseif(islogical(varargin{a}) && isscalar(varargin{a}))
        % Batch mode, don't display anything
        batch = varargin{a};
    elseif(ischar(varargin{a}) && any(strcmpi(varargin{a},{'nvt','nlx'})))
        % This is for a neuralynx file I believe
        nlx = true;
    else error('NEX2CMB:UnrecognizedInput', 'Argument %d unrecognized.',a+1);
    end
end


% If no input file name is provided then pass an empty string to next level
if(nargin < 1 || isempty(nexfile)); nexfile = ''; end


% Read the NEX file into a structure, or if it is already a structure,
% create the usable struct
if(isstruct(nexfile));
    nex = nexfile;
    if(isfield(nex,'nexfileinfo') && numel(nex.nexfileinfo)>=1 && ...
            isfield(nex.nexfileinfo,'name'))
        if(numel(nex.nexfileinfo)>1)
            names = char({nex.nexfileinfo(:).name}');
            samechar = all(diff(names+0,1,1)==0,1);
            nexfile = names(1,samechar);
            if(isempty(nexfile)); nexfile = nex.nexfileinfo(1).name; end
        else
            nexfile = nex.nexfileinfo.name;
        end
        if(~batch); fprintf(1,'Found NEX file name: %s\n',nexfile); end
    else
        nexfile = inputname(1);
        if(~batch); fprintf(1,'Using NEX file name: %s\n',nexfile); end
    end
else [nex,fid] = readNexFileM(nexfile);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  to interpret the video data %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine if and where the video data is stored in the NEX file


if(isfield(nex,'coords'))
    t = nex.coords(:,1);
    x = nex.coords(:,2);
    y = nex.coords(:,3);
    if(size(nex.coords,2)>=4)
        hd = nex.coords(:,[1 4]);
    end
    coords = [t x y];
    
elseif(isfield(nex,'rawLEDs'))
    if(~batch); fprintf(1,'Using ''rawLEDs'' instead of strobbed coordinates.\n'); end
    coords = nex.rawLEDs;
    
elseif(nlx)
    FilterSpec = {'*OrdCor.nvt', 'Neuralynx Video File (*.nvt)';
        '*.nvt', 'Neuralynx Video File (*.nvt)';
        '*', 'All Files'};
    [f, p] = uigetfile(FilterSpec,'Select a Neuralynx video (.nvt) file');
    if(f == 0); viddata = []; return; end
    f = fullfile(p,f);
    
  %  import CMBHOME.NL2Mat.Nlx2MatVT
    FieldSelect = [1 1 1 1 0 0];
    ExtractHeader = 1;
    ExtractMode = 1;
    [t, x, y, hd, head] = Nlx2MatVT(f, FieldSelect,ExtractHeader,ExtractMode);
    openedline = find(cellfun(@(x) ~isempty(x),strfind(head,'Time Opened')),1);
    mdyhms = regexp(head{openedline},...
        ': ([0-9]+)/([0-9]+)/([0-9]{4}) *At Time: ([0-9]+):([0-9]+):([0-9]+)','Tokens','Once');
    [mon,day,year,hour,min,sec] = mdyhms{:};
    nvtfiledate = datestr([year '-' mon '-' day ' ' hour ':' min ':' sec],31);
    
    if(~batch);
        fprintf(1,'Neuralynx Video Tracking File Date: %s\n',nvtfiledate);
    end
    
    t = t'/10^6; x = x'; y = y'; hd = [t hd'];
    coords = [t x y];

    
% this is likely the case we will be using.  This case gives you the option
% to import either nex or dvt coordinates.
elseif isfield(nex,'markers')
    coords = readNexCoords(nex);
    nex.coords = coords;
end


% THIS IS WHERE YOU WOULD MODIFY OR EDIT YOUR LED DATA
% fixnexxy or something to that degree


% if you have enough coordinates, you can import head direation and mean
% position data
if size(coords,2) <= 5

    % Find the mean position of all the LEDs
    if(exist('t','var')~=1 || exist('x','var')~=1 || exist('y','var')~=1)
    [t,x,y] = meanCoord(coords);
    end
    % Find the head direction based on the LEDs
    if(exist('hd','var')~=1 && size(coords,2)>=5)
        hd = headdirectionLED(coords);
        hd(:,2) = hd(:,2).*(180/pi);
    %     [hd, ~, LEDloc] = headdirectionLED(coords);
    else
    hd = nan(length(t),2);
    end
end


if exist('hd','var')
    if(max(abs(hd(:,2)))<=2*pi)
        warning('NEX2CMB:RadiansNotDegrees',...
        ['Head direction angles range from %f to %f, which appear to be radians.\n',...
        'However head direction should be given in degrees.'],...
        min(hd(:,2)),max(hd(:,2)));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% TO FIND AND ENTER EVENTS%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find out how many events we have.
if(isfield(nex,'events')); 
    nevin = size(nex.events,1);
else
    nevin = 0;
end

% Copy all events from the NEX file into the format needed by CMBHOME.
event = cell(0,2);
evcount = 0;
for ii = 1:nevin
    if(isempty(eventstokeep) || ismember(lower(nex.events{ii}.name),eventstokeep))
        newev = numel(nex.events{ii}.timestamps);
        event(evcount+1:evcount+newev,1) = {nex.events{ii}.name};
        event(evcount+1:evcount+newev,2) = num2cell(nex.events{ii}.timestamps);
        evcount = evcount+newev;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% NOW TIME TO MAKE THE OBJECT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% this creates the outfile for the cmb 
if isempty(nexfile)
    outfile = [fid(1:end-3) 'CMB.mat'];
else
    outfile = [nexfile(1:end-3) 'CMB.mat']; % why?????
end


% Create CMBHOME object.
viddata = CMBHOME.Session('name', outfile, 'b_ts', t,...
    'b_x',x, 'b_y', y, 'b_headdir', hd(:,2),...
    'event', event, 'raw_pos', 1, 'raw_headdir', 1, 'date_created', now, ...
    'epoch', t([1 end]), 'fs_video', 30, 'path_raw_data', nexfile);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% IF YOU WANT TO FIX THE POSITION DATA HERE%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% not sure if this will fix the head direction data though
%data = data.FixPos;

if(isempty(viddata.epoch)) 
    viddata.epoch = t([1 end]);
end
if(isfield(nex,'pxpercm')); 
    viddata.spatial_scale = 1/nex.pxpercm;
else
    viddata.spatial_scale = 0.25;
end
end

