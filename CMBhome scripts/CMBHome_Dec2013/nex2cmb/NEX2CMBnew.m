function [data, nex] = NEX2CMBnew(nexfile,varargin)

eventstokeep = {};
batch = false;
nlx = false;

% not necessary, this is used if you want to select events, or use a
% neuralynx nex file
% for a = 1:numel(varargin)
%     if(iscellstr(varargin{a}))
%         % 'eventstokeep' dictates which events from the NEX file are copied into
%         % the CMBHOME object. If empty, all events are kept. To exclude all events,
%         % include one element in 'eventstokeep' that won't match any other event,
%         % for example: eventstokeep = {''};
%         eventstokeep = lower(varargin{a});
%     elseif(islogical(varargin{a}) && isscalar(varargin{a}))
%         % Batch mode, don't display anything
%         batch = varargin{a};
%     elseif(ischar(varargin{a}) && any(strcmpi(varargin{a},{'nvt','nlx'})))
%         % This is for a neuralynx file I believe
%         nlx = true;
%     else error('NEX2CMB:UnrecognizedInput', 'Argument %d unrecognized.',a+1);
%     end
% end

% If no input file name is provided then pass an empty string to next level
if(nargin < 1 || isempty(nexfile)); nexfile = ''; end

% Read the NEX file into a structure.
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
else nex = readNexFileM(nexfile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  to interpret the video data %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find the position data stored in the NEX file.
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
    if(f == 0); data = []; return; end
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
elseif isfield(nex,'markers')
    coords = readNexCoords(nex);
    nex.coords = coords;
end

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
%vel = cat(1, 0, sqrt(diff(x).^2+diff(y).^2)) / mode(diff(t));
% Unless the LEDs were mounted on the cup, the HD is inaccurate.
% if(strcmp(LEDloc,'headstage')); hd(:,2) = nan(size(hd,1),1); end

% Create an output file name by strapping '.CMB.mat' to the end of the
% existing file name.
outfile = [nexfile(1:end-3) 'CMB.mat']; % why?????

% Find out how many events we have.
if(isfield(nex,'events')); nevin = size(nex.events,1);
else nevin = 0;
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

if exist('hd','var')
    if(max(abs(hd(:,2)))<=2*pi)
        warning('NEX2CMB:RadiansNotDegrees',...
        ['Head direction angles range from %f to %f, which appear to be radians.\n',...
        'However head direction should be given in degrees.'],...
        min(hd(:,2)),max(hd(:,2)));
    end
end

% Create CMBHOME object.
data = CMBHOME.Session('name', outfile, 'b_ts', t,...
    'b_x',x, 'b_y', y, 'b_headdir', hd(:,2),...
    'event', event, 'raw_pos', 1, 'raw_headdir', 1, 'date_created', now, ...
    'epoch', t([1 end]), 'fs_video', 30, 'path_raw_data', nexfile);
data = data.FixPos;

if(isempty(data.epoch)); data.epoch = t([1 end]); end
if(isfield(nex,'pxpercm')); data.spatial_scale = 1/nex.pxpercm;
else data.spatial_scale = 0.25;
end

% Get names for every neuron.
if isfield(nex,'neurons')
    names = cellfun(@(x) {x.name}, nex.neurons);
    
    % Get the tetrode numbers for each neuron.
    % This assumes that the tetrode number is the first block of consecutive
    % numbers within the tetrode name.
    tetnum = (str2double(regexp(names,'([0-9])+','once','match'))-1)/4;
    
    % Get the unit number for each neuron.
    % This assumes that the units are named such that the last letter of the
    % name denotes the unit number: a = 1, b = 2, c = 3, etc. This is how
    % Offline Sorter names units when using the 'Export to NEX' option.
    unitnum = cellfun(@(x) x(end)-96,names);
    
    
    %cineplex studio does not give tetrodes so delete repeats
    units = unique([tetnum unitnum],'rows');
    
    for jj = 1:size(units,1)
       keep(jj) = find(ismember( [tetnum unitnum],units(jj,:),'rows'),1,'first');
        
    end
    
    nex.neurons = nex.neurons(keep);
    tetnum = tetnum(keep);
    unitnum = unitnum(keep);
    for ii = 1:size(nex.neurons,1)
        if(~batch);
            fprintf(1,'%s is tetrode %d unit %d\n',names{ii},tetnum(ii),unitnum(ii)); 
        end
        spk_ts = nex.neurons{ii}.timestamps;
        spk_i = time2ind(t,spk_ts);
        spk_ts = spk_ts(spk_i>0);
        spk_i = spk_i(spk_i>0);
        
        data.spike(ii) = CMBHOME.Spike('i',spk_i,'ts',spk_ts,'tet_name',names{ii});
    end
end
end
