% This is a standard data format class for in vivo spiking, video tracking,
% and lfp data at the CMB at Boston University. See
% http://conte.bu.edu/repos/abogaard/+CMBHOME/README_CMBHOME.pdf and
% TutorialData.zip for a tutorial. 
%
% ACCESSING DOCUMENTATION
%
%   >> doc CMBHOME.Session; % for documentation on root objects
%   >> doc CMBHOME.LFP; % for documentation on LFP objects
%   >> doc CMBHOME.Spike; % for documentation on Spike objects (units)
%
% By default epochs work like ts>=epoch(1) & ts<=epoch(2). The hidden field
% root.ops = {'ge', 'le'} may change this
%
% All vectors are Nx1
%
% Version History
%
% v1 15 april 2010      presented and released
% v1.1 26 may 2010      changed root.spk_*(tetrode_i, cell_i) to
%                       root.spk_*(cel) to match convention
% v1.2 18 june 2010     added epoch_group to indicate epochs which indicate
%                       the same type of event
% v1.3 8 july 2010      added notes property
% v1.4 18 august 2010   added active LFP file property and how we handle
%                       LFP objects
% v1.5 11 jan 2011      added spk property, as well as user-defined "cel"
%                       property, updated behavior of root.vel. (if
%                       root.b_vel exists, it uses this. if not, simple
%                       derivative is used). also added the
%                       root.name_formatted field so that now the objects
%                       work in a relative naming convention if desired
% v2   14 oct 2011      greatly improved speed and inter-dependencies
%                       between properties
% andrew robert bogaard 2/10-9/11

classdef Session
    
    properties (Hidden)
        
        ops = {'ge', 'le'}  % operators for epochs behavior
        
    end
    
    properties (Access=private)
       
        p_ind               % same as root.ind, but set by set.epoch
        p_lfp_ind
        p_cel_ind           % indices in root.b_ts and root.b_lfp(active_lfp).ts for each cell spike. set with set.cel
        p_cel_spkind        % indices in root.spike(cel(1), cel(2)).ts for each spike in the epoch. set with set.cel
        p_cel_lfp_ind
        
    end
    
    properties
        
        date_created        % date and time created as a serial date number (see 'doc now')
        notes               % cell array of notes and dates {datestamp, str_note}
        spike               % object array of class 'Spike'  
        name                % 1xM cellstr of directory and filename of current object. see name_formatted and Save. 
                            % If you want to use a relative path, use a cell array of strings like 
                            % {'dir1', 'dir1', 'filename.mat'}
        cell_thresh = [0,0] % threshold of minimum spike rate [minimum spike frequency, check drift]
        fs_video            % sampling rate from video (samples/sec)
        raw_pos=1           % boolean, is this uncorrected position?
        raw_headdir=1       % boolean, is this uncorrected head direction?
        event = cell(1, 2)  % event flags, Nx2 cell array, {'label', time}            
        path_lfp            % cell array of fnames to lfp recording data
        path_raw_data       % path to original recording data
        spatial_scale=.5    % cm/pixel
        
        epoch_group         % vector, Nx1, N is number of epochs with integer indicating like-epochs
        user_def            % property that may be whatever user defines, ex. struct of relevant
        version = 2         % version control
        b_lfp               % object array of structs with fields ts, signal, path_lfp_index ex
                   
    end
    
    properties (Access=private, Hidden)
        
        p_epoch
        p_cel
        p_active_lfp
        p_b_x               
        p_b_y                 
        p_b_vel               
        p_b_headdir           
        p_b_ts               
        p_b_myvar
        
    end
    
    properties (Dependent)
        
        epoch               % array, Nx2, start and stop times that exists (can be inf or -inf)
        cel                 % Nx2 array of [tetrode, cell index]
        active_lfp          % index to root.path_lfp indicating which LFP is loaded (or to be loaded)
        
    end
    
    properties (Dependent, Hidden)
    
        b_x                 % x position
        b_y                 % y position
        b_vel               % user defined velocity, if undefined, the derivative of b_x and b_y is used
        b_headdir           % head direction
        b_ts                % timestamps
        b_myvar
    
    end
    
    properties (Dependent=true, SetAccess=private)
        
        vel                 % velocity, within start and stop epoch (pixels/sec)
        x                   % x position (pixels)
        y                   % y position (pixels)
        ind                 % index in root.b_ts
        ts                  % timestamps (seconds)
        myvar
        b_epoch_group       % vector, Nx1, N is number of epochs with integer indicating like-epochs
        headdir             % head direction, (degrees, -180:180)
        cells               % Nx2 array of cells which meet specs_cell above,  [tetrode, cell;...]
        active_event        % an Nx2 element cell array with event labels (if exist) that match current epoch times.
                            % N is the number of rows in root.epochs
        lfp                 % object array of LFP objects of fields ts and signal
        spk                 % struct with fields representing data at spike times of root.cel
        name_formatted      % this is the root.name field formatted for the current platform
        cel_x
        cel_y
        cel_vel
        cel_ts
        cel_headdir
        cel_i
        cel_theta
        cel_thetaphase
        cel_myvar
                
    end
    
    methods (Static)
        
        [t,x,y,vx,vy,ax,ay] = KalmanVel(posx,posy,post,order,Q,R);
        
        function root = MergeSessions(cellarray_filename)
        % (1) root = CMBHOME.Session.MergeSessions
        % (2) root = CMBHOME.Session.MergeSessions(cellarray_filenames)
        %
        % Merges multiple Session objects by stitching together time
        % indeces. Saves original Session times as events
        %
        % (1) Prompts user to select files to merge, stiches them together
        % in the order selected
        % (2) Accepts a cell array of file names to be merged. Stitches
        % them together in the order of cellarray_filenames
        %
        % andrew 24 may 2010
        
        warning('CMBH:notify', '%s', 'This function is not tested, and could be tweaked so that it is more useful. Come see me if you have comments or suggestions, andrew.');
        
        import CMBHOME.*
        f_ind=1;
        if ~exist('cellarray_filenames', 'var')
            while 1 % loop until user cancels, and we get no files

                [load_files, base_path] = uigetfile('*.mat','Please select Session object mat files to merge. Exit to finish.',pwd, 'MultiSelect', 'on');

                if iscell(load_files)
                        for i = 1:length(load_files)
                            fopen(f_ind) = {fullfile(base_path, load_files{i})};
                            f_ind = f_ind+1;
                        end
                elseif ischar(load_files)
                        fopen(f_ind) = {fullfile(base_path, load_files)};
                        f_ind = f_ind+1;
                else
                    break; % exit loop collecting files
                end

            end
        else
            fopen = cellarray_filenames;
        end
        
        t_start = 0;
        
        ind_offset = 0;
        
        root = Session('path_lfp', 'Merged Sessions, see root.path_raw_data for original file names', ...
            'path_raw_data', fopen, 'date_created', now, 'name', 'Merged Sessions, see root.path_raw_data.');
        
        for i = 1:length(fopen)
            
            tmp = load(fopen{i});
            
            offset = tmp.root.b_ts(1) - t_start;
            
            if i==1, 
                root.fs_video = tmp.root.fs_video; 
                root.spatial_scale = tmp.root.spatial_scale;
            end
            
            root.b_ts = cat(1, root.b_ts, tmp.root.b_ts-offset);
            
            tmp.root.event(:,2) = num2cell([tmp.root.event{:,2}]-offset);
            
            root.b_x = cat(1, root.b_x, tmp.root.b_x);
            root.b_y = cat(1, root.b_y, tmp.root.b_y);
            root.b_headdir = cat(1, root.b_headdir, tmp.root.b_headdir);
            root.event = cat(1, root.event, tmp.root.event, {[fopen{i} ' start'], t_start}, {[fopen{i} ' end'], t_start + (tmp.root.b_ts(end)-tmp.root.b_ts(1))});
            
            tmp.root.cell_thresh = [ 0 0 ];
            
            for j = 1:size(tmp.root.cells,1)
                
                tmp.root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).ts = tmp.root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).ts - offset; % shift indices and timestamps
                tmp.root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).i = tmp.root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).i + ind_offset;
                
                if all(size(root.spike)>=tmp.root.cells(j,:))
                    root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).ts = cat(1, root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).ts, tmp.root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).ts);
                    root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).ts = cat(1, root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).i, tmp.root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)).i);
                else
                    root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2)) = tmp.root.spike(tmp.root.cells(j,1), tmp.root.cells(j,2));
                end
                
            end
            
            ind_offset = ind_offset + length(tmp.root.b_ts)+1;
            t_start = t_start + (tmp.root.b_ts(end)-tmp.root.b_ts(1))+1; % 1 second buffer
            
        end
        
        end
    end
    
    methods
                
        function self = Session(varargin)
            
            import CMBHOME.Spike
            
            p = inputParser;

            p.addParamValue('name',          'default session name', @(x) ischar(x));
            p.addParamValue('b_x',          [], @(x) any(size(x)<=1)); 
            p.addParamValue('b_y',          [], @(x) any(size(x)<=1)); 
            p.addParamValue('b_headdir',        [], @(x) any(size(x)<=1)); 
            p.addParamValue('b_ts',         [], @(x) any(size(x)<=1)); 
            p.addParamValue('b_lfp',         CMBHOME.LFP, @(x) isstruct(x));
            p.addParamValue('fs_video',     [], @(x) length(x)==1);
            p.addParamValue('raw_pos',      [], @(x) length(x)==1);
            p.addParamValue('raw_headdir',      [], @(x) length(x)==1);
            p.addParamValue('event',       cell(1,2), @(x) iscell(x));
            p.addParamValue('rem_ts',       [], @(x) length(x)==1);
            p.addParamValue('ripple_ts',    [], @(x) length(x)==1);
            p.addParamValue('path_lfp',     []);
            p.addParamValue('path_raw_data',[]);
            p.addParamValue('date_created', now, @(x) isnumeric(x));
            p.addParamValue('spatial_scale',.5, @(x) length(x)==1);
            p.addParamValue('epoch',        [0 1], @(x) numel(x)==2);
            p.addParamValue('spike',        Spike, @(x) strcmp(class(x), 'CMBHOME.Spike'));
                 
            p.parse(varargin{:});
            
            self.name = p.Results.name;
            self.notes = { date, ['Session object created, and saved to ', p.Results.name] };
            self.b_x = p.Results.b_x;
            self.b_y = p.Results.b_y;
            self.b_headdir = p.Results.b_headdir;
            self.b_ts = p.Results.b_ts;
            %self.b_lfp = p.Results.b_lfp;
            self.fs_video = p.Results.fs_video;
            %self.fs = p.Results.fs;
            self.raw_pos = p.Results.raw_pos;
            self.raw_headdir = p.Results.raw_headdir;
            self.event = p.Results.event;
            self.path_lfp = p.Results.path_lfp;
            self.path_raw_data = p.Results.path_raw_data;
            self.date_created = p.Results.date_created;
            self.spatial_scale = p.Results.spatial_scale;

            self.epoch = p.Results.epoch;
            self.spike = p.Results.spike;
        end   
        
        
        %% DEPENDENT STATE VARIABLES
        
        function self = set.epoch(self, epoch)
        % makes sure that all epochs are valid periods of time, and that
        % the array is formatted correctly. 
               
            if isempty(self.b_ts), error('root.b_ts must be set to use epochs'); end
        
            tmin = self.b_ts(1);
            tmax = self.b_ts(end);
            
            if numel(epoch)==2, epoch = epoch(:)'; end

            if size(epoch,2)~=2
                error('Epoch must be a matrix with columns [start times, stop times]');
            end

            if any(epoch(:,2)-epoch(:,1)<=0) % check that all starts and stops are start<stop
                error('epoch must be in format [tstart, tstop]');
            end
            
            epoch(epoch==inf) = self.b_ts(end);
            epoch(epoch==-inf) = self.b_ts(1);
            
%             epoch2 = CMBHOME.Utils.SelectEpochs(epoch, [tmin tmax], 1); % check that all epochs land within root.b_ts 
%             
%             if ~all(epoch2==epoch), disp('Warning, not all epochs were contained within root.b_ts. Those outside the recording time were removed, and those overreaching bounds were clipped'); end
%             
%             epoch = epoch2;


            if all(size(epoch)==size(self.p_epoch))
                
                if all(epoch==self.p_epoch) % you are reassigning the same thing!! dont run all the scripts below
                    
                    return
                    
                end
                
            end

            self.p_epoch = epoch;
            
            self.p_ind = IsolateEpoch(self.b_ts, [1:length(self.b_ts)]', self.epoch);
            
            self = SetLFPInds(self); % update the aligned lfp inds
            
            self = SetCelInds(self); % update the indices for the cells which are set
            
        end  
        
        function epoch = get.epoch(self)
            
            epoch = self.p_epoch;
            
        end
        
        function self = set.active_lfp(self, active_lfp)
            
            if active_lfp==self.p_active_lfp; return; end % its already set, don't rerun the indexing
            
            if active_lfp>0 & numel(active_lfp) == 1 & active_lfp<=length(self.b_lfp)
                self.p_active_lfp = active_lfp;
            elseif isempty(active_lfp)
                self.p_active_lfp = [];
            else
                disp('root.active_lfp not set because it was not a valid index for root.b_lfp')
                self.p_active_lfp = [];
            end
            
            self = SetCelInds(self);
            
        end
        
        function active_lfp = get.active_lfp(self)
            
            active_lfp = self.p_active_lfp;
            
        end
        
        function self = set.cel(self, cel)
        % sets the root.cel property (Nx2, [tetrode, cell_ind])

            cells = ValidCells(self);
        
            if isempty(cel), self.p_cel = []; return, end
        
            if size(cel, 2)~=2, warning('CMBH:error', 'root.cel must be Nx2, [tetrode, cell_ind]'); return; end
        
            if ~isnumeric(cel), warning('CMBH:error', 'root.cel must be Nx2, [tetrode, cell_ind]'); return; end
            
            if any(~ismember(cel, cells, 'rows'))
                warning('CMBH:error', 'error, some of root.cel are not cells, or do not pass root.cell_thresh. removed the following cells:');
                disp(num2str(cel(find(~ismember(cel, cells, 'rows')),:)));
                cel(find(~ismember(cel, cells, 'rows')),:) = [];
            end
            
            if all(size(cel)==size(self.p_cel))
                
                if all(cel==self.p_cel) % you're trying to reset the same cells, so dont rerun the indexing
                    
                    return
                    
                end
                
            end
            
            self.p_cel = cel;
            
            self = SetCelInds(self); % update the indices for the cells which are set
            
        end
        
        function cel = get.cel(self)
            
            cel = self.p_cel;
            
        end
        
        
        %% DEPENDENT BASE VARIABLES
        
        function self = set.b_x(self, b_x)
        % sets p_b_x. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_x) && ~isvector(b_x)) & ~isempty(b_x)
                error('X data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_x)~=l) & ~isempty(b_x)
                
                error('b_x must be the same length as all other base variables');
                
            end
            
            self.p_b_x = b_x(:);
            
        end
        
        function self = set.b_y(self, b_y)
        % sets p_b_y. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_y) && ~isvector(b_y)) & ~isempty(b_y)
                error('y data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_y)~=l) & ~isempty(b_y)
                
                error('b_y must be the same length as all other base variables');
                
            end
            
            self.p_b_y = b_y(:);
            
        end
        
        function self = set.b_vel(self, b_vel)
        % sets p_b_vel. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_vel) && ~isvector(b_vel)) & ~isempty(b_vel)
                error('vel data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_vel)~=l) & ~isempty(b_vel)
                
                error('b_vel must be the same length as all other base variables');
                
            end
            
            self.p_b_vel = b_vel(:);
            
        end
        
        function self = set.b_headdir(self, b_headdir)
        % sets p_b_headdir. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_headdir) && ~isvector(b_headdir)) & ~isempty(b_headdir)
                error('headdir data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_headdir)~=l) & ~isempty(b_headdir)
                
                error('b_headdir must be the same length as all other base variables');
                
            end
            
            self.p_b_headdir = b_headdir(:);
            
        end
        
        function self = set.b_ts(self, b_ts)
        % sets p_b_ts. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_ts) && ~isvector(b_ts)) & ~isempty(b_ts)
                error('ts data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_ts)~=l) & ~isempty(b_ts)
                
                error('b_ts must be the same length as all other base variables');
                
            end
            
            self.p_b_ts = b_ts(:);
            
        end
        
        function self = set.b_myvar(self, b_myvar)
        % sets p_b_myvar. must be either empty or the same length as other base variables
        
            if (~isnumeric(b_myvar) && ~isvector(b_myvar)) & ~isempty(b_myvar)
                error('y data must be numerical vector');
            end
            
            l = CheckBaseVarLength(self);
            
            if (l~=0 && length(b_myvar)~=l) & ~isempty(b_myvar)
                
                error('b_myvar must be the same length as all other base variables');
                
            end
            
            self.p_b_myvar = b_myvar(:);
            
        end
        
        function b_x = get.b_x(self)
            
            b_x = self.p_b_x;
            
        end
        
        function b_y = get.b_y(self)
            
            b_y = self.p_b_y;
            
        end
        
        function b_ts = get.b_ts(self)
            
            b_ts = self.p_b_ts;
            
        end
        
        function b_vel = get.b_vel(self)
            
            b_vel = self.p_b_vel;
            
        end
        
        function b_headdir = get.b_headdir(self)
            
            b_headdir = self.p_b_headdir;
            
        end
        
        function b_myvar = get.b_myvar(self)
            
            b_myvar = self.p_b_myvar;
            
        end
        
        
        %% DEPENDENT NON-BASE VARIABLES
        
        function cells = get.cells(self)
            cells = ValidCells(self);
        end
        
        function name_formatted = get.name_formatted(self)
        % returns the filename of the current object formatted for the current platform
        
            if isempty(self.name)

                name_formatted = 'CMBHno_name.mat';

            elseif iscellstr(self.name)

                name_formatted = fullfile(self.name{:});

            else

                name_formatted = self.name;

            end
            
        end
        
        function active_event = get.active_event(self)
        % e = root.active_event;
        %
        % Returns strings which correspond to the start and stop labels for
        % each epoch time.
        %
        % RETURNS
        %
        % e -> Nx2 cell array of strings (same size as root.epoch, or if
        % root.epoch_group is set, N is the number of groups, and it 
        % defaults to the group's first label)
        %
        % andrew 18 june 2010
        
            if size(self.event,2)<2
                
                active_event = {'No Events', 'No Events'};
            
            else
            
                active_event = cell(size(self.epoch,1), 2);
                counter = 1;
                for i = self.b_epoch_group'
                    
                    inds = find(self.b_epoch_group==i);
                    
                    if sum(vertcat(self.event{:,2})==self.epoch(inds(1), 1))==1
                        str1 = self.event{vertcat(self.event{:,2})==self.epoch(inds(1), 1),1};
                    else
                        str1 = ['Unknown t = ' num2str(self.epoch(inds(1), 1))];
                    end

                    if sum(vertcat(self.event{:,2})==self.epoch(inds(1), 2))==1
                        str2 = self.event{vertcat(self.event{:,2})==self.epoch(inds(1), 2),1};
                    else
                        str2 = ['Unknown t = ' num2str(self.epoch(inds(1), 2))];
                    end

                    active_event(counter, 1) = {str1};
                    active_event(counter, 2) = {str2};
                    
                    counter = counter+1;
                end
            
            end
            
        end
        
        function vel = get.vel(self)
        % returns pixels/sec
        
            if ~isempty(self.b_vel) % if there is a user defined velocity
                
                vel = [];

                if iscell(self.p_ind)
                    vel = cellfun(@(c) self.b_vel(c), self.p_ind, 'UniformOutput', false);
                elseif ~isempty(self.b_vel)
                    vel=self.b_vel(self.p_ind);
                end
                
            else
            
                if ~isempty(self.x) & length(self.x)==length(self.y)

                    if iscell(self.x)

                        vel = cell(length(self.x),1);

                        for i = 1:length(self.x)
                            vel{i} = cat(1, 0, sqrt(diff(self.x{i}).^2+diff(self.y{i}).^2))*self.fs_video;
                        end

                    else
                        vel = cat(1, 0, sqrt(diff(self.x).^2+diff(self.y).^2))*self.fs_video;
                    end

                else

                    warning('CMBH:error', 'Incoherent x and y tracking, or no tracking.');

                    vel = [];

                end
            end
        end
                
        function ind = get.ind(self)
        % These are indices in the base vectors (b_*) used to return data dynamic fields (ex. root.x = root.b_x(root.ind))
        %
        % This is useful if you have a user defined vector that isnt
        % included in the set of properties already defined. For example, a
        % linearized position in the TMaze could be saved in
        % root.user_def.linear_t. Then for any epoch the linear position on
        % the track is root.user_def.linear_t(root.ind), so long as
        % length(linear_t)==length(root.b_ts)        
        
            ind = self.p_ind;
                        
        end
        
        function x = get.x(self)
        % x position in pixels
            
            x = [];
                        
            if iscell(self.p_ind)
                x = cellfun(@(c) self.b_x(c), self.p_ind, 'UniformOutput', false);
            elseif ~isempty(self.b_x)
                x=self.b_x(self.p_ind);
            else
                warning('CMBH:notify', '%s', 'No x data for this session.');
            end
            
        end
        
        function y = get.y(self)
        % y position in pixels

            y = [];
            
            if iscell(self.p_ind)
                y = cellfun(@(c) self.b_y(c), self.p_ind, 'UniformOutput', false);
            elseif ~isempty(self.b_y)
                y=self.b_y(self.p_ind);
            else
                warning('CMBH:notify', '%s', 'No y data for this session');
            end
            
        end
        
        function headdir = get.headdir(self)
        % headdirection in degrees (although the user could input radians, 
        % all Toolbox functionality assumes degrees)  

            headdir = [];
                     
                if iscell(self.p_ind)
                    headdir = cellfun(@(c) self.b_headdir(c), self.p_ind, 'UniformOutput', false);
                elseif ~isempty(self.b_headdir)
                    headdir=self.b_headdir(self.p_ind);
                else
                    warning('CMBH:notify', '%s', 'Head Direction vector does not exist for this session');
                end
        end
        
        function ts = get.ts(self)
        % timestamps in seconds    
            
            if iscell(self.p_ind)
                ts = cellfun(@(c) self.b_ts(c), self.p_ind, 'UniformOutput', false);
            else
                ts=self.b_ts(self.p_ind);
            end
            
        end
        
        function myvar = get.myvar(self)
        % epoch data from b_myvar

            if isempty(self.b_myvar)
                myvar = [];
                return
            elseif length(self.b_myvar)~=length(self.b_ts) && ~isempty(self.b_myvar)
                disp('root.b_myvar must be the same length as root.b_ts');
                myvar = [];
                return

            end

            myvar = [];

            if iscell(self.p_ind)
                myvar = cellfun(@(c) self.b_myvar(c), self.p_ind, 'UniformOutput', false);
            elseif ~isempty(self.b_myvar)
                myvar=self.b_myvar(self.p_ind);
            else
                warning('CMBH:notify', '%s', 'No myvar data for this session');
            end

        end
        
        function lfp = get.lfp(self)
            % check to see if b_lfp is populated. if not, check that there
            % is a file by the index at root.active_lfp
            
            % if there is, load lfp it
            
            % then call self.lfp.signal and self.lfp.time and get chucks as
            % per self.epoch
            
            if exist('self.b_lfp') % check if it's loaded yet
                if length(self.b_lfp)<self.active_lfp
                    
                    error('root.active_lfp exceeds length of root.b_lfp');
                    
                end
            elseif ~isempty(self.active_lfp)
                
                signal = self.b_lfp(self.active_lfp).signal;
                ts = self.b_lfp(self.active_lfp).ts;
                              
                b_theta = self.b_lfp(self.active_lfp).b_theta;

                b_theta_phase = self.b_lfp(self.active_lfp).b_theta_phase;

                b_theta_amplitude = self.b_lfp(self.active_lfp).b_theta_amplitude;
                
                if isempty(b_theta), b_theta = nan(length(signal), 1); end
                
                if isempty(b_theta_phase), b_theta_phase = nan(length(signal), 1); end
                
                if isempty(b_theta_amplitude), b_theta_amplitude = nan(length(signal), 1); end
            
            else
                
                lfp = [];
                return
            
            end
            %% grab data from epochs
            
            if size(self.epoch,1)>1
                
                tmp_signal = cell(size(self.epoch,1), 1);
                tmp_ts = cell(size(self.epoch,1), 1);
                tmp_b_theta = cell(size(self.epoch,1), 1);
                tmp_b_theta_phase = cell(size(self.epoch,1), 1);
                tmp_b_theta_amplitude = cell(size(self.epoch,1), 1);
            
                for i = 1:size(self.epoch,1)
                    
                    inds = ts>=self.epoch(i,1) & ts<=self.epoch(i,2);
                                        
                    tmp_signal{i} = signal(inds);
                    tmp_ts{i} = ts(inds);
                    tmp_b_theta{i} = b_theta(inds);
                    tmp_b_theta_phase{i} = b_theta_phase(inds);
                    tmp_b_theta_amplitude{i} = b_theta_amplitude(inds);
                    
                end
                
                signal = tmp_signal;
                ts = tmp_ts;
                b_theta = tmp_b_theta;
                b_theta_phase = tmp_b_theta_phase;
                b_theta_amplitude = tmp_b_theta_amplitude;
                
            else
                
                inds = ts>=self.epoch(1) & ts<=self.epoch(2);
                
                signal = signal(inds);
                ts = ts(inds);
                b_theta = b_theta(inds);
                b_theta_phase = b_theta_phase(inds);
                b_theta_amplitude = b_theta_amplitude(inds);
                
            end
            
            user_def = self.b_lfp(self.active_lfp).user_def;
            
            fs = self.b_lfp(self.active_lfp).fs;
            
            channel_name = self.b_lfp(self.active_lfp).channel_name;
            
            lfp = CMBHOME.LFP(signal, ts, fs, channel_name, b_theta, b_theta_phase, b_theta_amplitude, user_def);
            
        end
        
        function b_epoch_group = get.b_epoch_group(self)
        % Checks to see if root.epoch_group is set by user. If it is, and 
        % is valid (length=number of epochs), then root.b_epoch_groups ==
        % root.epoch_groups
        
            if length(self.epoch_group) == size(self.epoch,1) && all(self.epoch_group>0)

                b_epoch_group = self.epoch_group;

            else

                b_epoch_group = 1:size(self.epoch,1);

                b_epoch_group = b_epoch_group';

            end            
            
        end
        
        function spk = get.spk(self)
        % builds the dependent property root.spk, which is a struct with fields for 
        % all relevant data that occurred at spike times for root.cel
        
            disp('sadly, the root.spk. syntax is terribly slow. check out root.cel_x etc. seriously. rewrite this part of your code (sorry).');
        
            if isempty(self.cel), spk = []; return; end
            
            if ~isempty(self.b_vel), vel = self.b_vel; 
            else vel = cat(1, 0, sqrt(diff(self.b_x).^2+diff(self.b_y).^2))/self.fs_video^-1;
            end
                          
            if isempty(self.b_myvar), self.b_myvar = nan(length(self.b_ts),1); end

            ind_lfp = self.active_lfp;
        
            cel = self.cel;

            if size(self.epoch,1)>1 || size(cel,1)>1
                spk.x = cell(size(self.epoch,1) , size(cel,1));
                spk.y = cell(size(self.epoch,1) , size(cel,1));
                spk.myvar = cell(size(self.epoch,1) , size(cel,1));
                spk.ts = cell(size(self.epoch,1) , size(cel,1));
                spk.vel = cell(size(self.epoch,1) , size(cel,1));
                spk.i = cell(size(self.epoch,1) , size(cel,1));
                spk.headdir = cell(size(self.epoch,1) , size(cel,1));
                spk.theta = cell(size(self.epoch,1) , size(cel,1));
            else
                spk.x = [];
                spk.y = [];
                spk.myvar = [];
                spk.ts = [];
                spk.vel = [];
                spk.i = [];
                spk.headdir = [];
                spk.theta = [];
            end
            
            for j = 1:size(cel,1)
            
                [i, i2] = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

                if iscell(i)
                    spk.x(:,j) = cellfun(@(c) self.b_x(c), i, 'UniformOutput', false);
                    spk.y(:,j) = cellfun(@(c) self.b_y(c), i, 'UniformOutput', false);
                    spk.vel(:,j) = cellfun(@(c) vel(c), i, 'UniformOutput', false);
                    spk.myvar(:,j) = cellfun(@(c) self.b_myvar(c), i, 'UniformOutput', false);
                    spk.headdir(:,j) = cellfun(@(c) self.b_headdir(c), i, 'UniformOutput', false);
                    spk.i(:,j) = i;
                    spk.ts(:,j) = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).ts(c), i2, 'UniformOutput', false);
                    
                    if ind_lfp
                        i2 = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).i_lfp{ind_lfp}(c), i2, 'Unif', 0);
                        spk.theta(:,j) = cellfun(@(c) self.b_lfp(ind_lfp).theta_phase(c), i2, 'UniformOutput', false);
                    end
                        
                elseif size(cel,1)>1
                    spk.x{1,j} = self.b_x(i);
                    spk.y{1,j} = self.b_y(i);
                    spk.vel{1,j} = vel(i);
                    spk.myvar{1,j} = self.b_myvar(i);
                    spk.headdir{1,j} = self.b_headdir(i);
                    spk.i{1,j} = i;
                    spk.ts{1,j} = self.spike(cel(j,1), cel(j,2)).ts(i2);
                    
                    if ind_lfp
                        i2 = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).i_lfp{ind_lfp}(c), i2, 'Unif', 0);
                        spk.theta{1,j} = self.b_lfp(ind_lfp).theta_phase(i2);
                    end
                    
                else 
                    spk.x = self.b_x(i);
                    spk.y = self.b_y(i);
                    spk.vel = vel(i);
                    spk.myvar = self.b_myvar(i);
                    spk.headdir = self.b_headdir(i);
                    spk.i = i;
                    spk.ts = self.spike(cel(j,1), cel(j,2)).ts(i2);
                    
                    if ind_lfp
                        i2 = self.spike(cel(j,1), cel(j,2)).i_lfp{ind_lfp}(i2);
                        spk.theta = self.b_lfp(ind_lfp).theta_phase(i2);
                    end
                end

            end
                    
        end
        
        function cel_x = get.cel_x(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_x = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            cel_x = cellfun(@(c) self.b_x(c), self.p_cel_ind, 'unif', 0);
            
        end
                        
        function cel_y = get.cel_y(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_y = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            cel_y = cellfun(@(c) self.b_y(c), self.p_cel_ind, 'unif', 0);
            
        end

        function cel_vel = get.cel_vel(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_vel = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            if ~isempty(self.b_vel),
            
                cel_vel = cellfun(@(c) self.b_vel(c), self.p_cel_ind, 'unif', 0);
            
            else
                
                vel = cat(1, NaN, sqrt(diff(self.b_x).^2+diff(self.b_y).^2))*self.fs_video;
                cel_vel = cellfun(@(c) vel(c), self.p_cel_ind, 'unif', 0);
                
            end

        end
            
        function cel_ts = get.cel_ts(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_ts = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            for i = 1:size(self.cel,1)
                
                cel = self.cel(i,:);
            
                cel_ts(:,i) = cellfun(@(c) self.spike(cel(1), cel(2)).ts(c), self.p_cel_spkind(:,i), 'unif', 0);
            
            end
                
        end

        function cel_headdir = get.cel_headdir(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_headdir = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            cel_headdir = cellfun(@(c) self.b_headdir(c), self.p_cel_ind, 'unif', 0);
            
        end

        function cel_i = get.cel_i(self)
            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_i = self.p_cel_ind;
                
                return
                
            end
            
            cel_i = self.p_cel_ind;
            
        end

        function cel_myvar = get.cel_myvar(self)
                            
            if isempty(self.p_cel_ind)
                
                disp('Set root.cel.'); 
            
                cel_myvar = cell(size(self.epoch,1), 1);
                
                return
                
            end
            
            cel_myvar = cellfun(@(c) self.b_myvar(c), self.p_cel_ind, 'unif', 0);
            
        end
            
        function cel_theta = get.cel_theta(self)
           
            if isempty(self.cel) % if no cel is set, return empty
                disp('Set root.cel');
                cel_theta = cell(size(self.epoch,1), 1);
                return
            end
            
            if isempty(self.active_lfp) % if no active lfp set 
                disp('Set root.active_lfp');
                cel_theta = cell(size(self.epoch,1), size(self.cel, 1));
                return                
            end
            
            cel_theta = cellfun(@(c) self.b_lfp(self.active_lfp).theta(c), self.p_cel_lfp_ind, 'unif', 0);
                
        end
        
        function cel_thetaphase = get.cel_thetaphase(self)
           
            if isempty(self.cel) % if no cel is set, return empty
                disp('Set root.cel');
                cel_thetaphase = cell(size(self.epoch,1), 1);
                return
            end
            
            if isempty(self.active_lfp) % if no active lfp set 
                disp('Set root.active_lfp');
                cel_thetaphase = cell(size(self.epoch,1), size(self.cel, 1));
                return                
            end
            
            cel_thetaphase = cellfun(@(c) self.b_lfp(self.active_lfp).theta_phase(c), self.p_cel_lfp_ind, 'unif', 0);
                
        end
                
        
        %% INDEPENDENT PROPERTIES
        
        function self = set.cell_thresh(self, cell_thresh)
            
            if isnumeric(cell_thresh)
                
                if numel(cell_thresh)==1
                    cell_thresh = [cell_thresh(1), 0];
                elseif numel(cell_thresh)==2
                    if cell_thresh(2)~=1 && cell_thresh(2)~=0
                        warning('CMBH:error', 'cell_thresh(2) must be boolean: average Fmin (0) or sustained Fmin (1)? Defaulting to 0.')
                        cell_thresh(2) = 0;
                    end
                else
                    warning('CMBH:error', 'Improper number of elements on cell_thresh');
                end
                
                switch cell_thresh(2)
                case 1
                    str_cell = ' Hz continuously throughout the entire session';
                case 0
                    str_cell = ' Hz on average during the session';
                end
                
                self.cell_thresh = cell_thresh;

                warning('CMBH:notify', '%s', ['Cell Threshold: F_min is ' num2str(self.cell_thresh(1)) str_cell '.']);
                
            end
        end
        
        function self = set.epoch_group(self, epoch_group)
           
            if all(epoch_group>0) && all(rem(epoch_group, 1)==0) && isvector(epoch_group)
               
                self.epoch_group = epoch_group(:);
                
            end
            
        end
                        
        function self = set.event(self, event)
            if ~iscell(event)
                error('session.event must be cell array')
            else
                self.event=event;
            end
            
            if size(event, 2)==2 && size(event,1)>1 % sort events if numeric
            
                isnumbers = cellfun(@isnumeric, event(:,2));
                isemptys = cellfun(@isempty, event(:,2));
                
                if any(~isnumbers) || any(isemptys)
                    
                    warning('CMBH:notify', '%s', 'There are values in the timestamps column (2) that are not numbers. They will be deleted');
                    
                    event(~isnumbers,:) = [];
                    event(isemptys,:) = [];
                    
                end
                
                ts = [event{:,2}];
                
                if isnumeric(ts)
                   
                    [trash, ind] = sort(ts); 
                    
                    self.event = event(ind,:);
                    
                end
                
            elseif isempty(event)
                
                self.event = {'Session start', self.b_ts(1); 'Session stop', self.b_ts(end) };
                
            end
            
        end
        
        function self = set.path_lfp(self, path_lfp)
            if isempty(path_lfp)
                self.path_lfp = [];
                return;
            end
            
            if ~iscell(path_lfp), path_lfp = {path_lfp}; end
            
%             if iscell(path_lfp)
%                 ext = cellfun(@(c) c(end-3:end), path_lfp, 'UniformOutput', false);
%             else
%                 ext = path_lfp(end-3:end);
%             end
            
            %if sum([ismember(ext,'.mat'), ismember(ext, '.plx'), ismember(ext, '.ncs')]) == length(path_lfp)
                self.path_lfp = path_lfp;
            %else
            %    error('Error seeting path_lfp: lfp files must be .plx, .ncs, or .mat files. We can update this if need be.');
            %end
              
        end
       
        function self = set.b_lfp(self, b_lfp)
           
            issue_warning = 0;
            
            if strcmp(class(b_lfp), 'CMBHOME.LFP')
            
                for i = 1:length(b_lfp)
                    
                    if any(unique(diff(b_lfp(i).ts)).^-1 < .9*b_lfp(i).fs), issue_warning = 1; end

                end
                
            end
            
            if issue_warning, warning('CMBH:notify', '%s', 'There are timestamps with uneven sampling intervals. Make sure this is expected for this session, and to use root.epoch properly.'); end
        
            self.b_lfp = b_lfp; 
            
        end
        
        
        %% METHODS
        
        function [spk_xcorr, lags, epoch, std_error] = spk_xcorr(self, cells, max_lag, t_bin, average_epochs, norm)
        % [cor, lag, epoch, std_error] = root.spk_xcorr(cel, max_lag,t_bin, average_epochs);
        %
        % Returns the unbiased crosscorrelation, or autocorrelation for cell(s) in
        % cel
        %
        % ARGUMENTS
        %   cells               a Nx2 (N=1 or 2) array of cells for which to
        %                       calculate xcorr. if N=1, than autocorrelation is returned
        %   max_lag             maximum lag for which xcorr is calculated. epochs that are not this long are removed
        %   t_bin               binsize in seconds
        %   average_epochs      1 or 0 (default 0). If 1, cor is an array
        %                       calculated by the average of correlations for each epoch
        %                       weighted by epoch duration. If 0, cor is a
        %                       cell array of correlation vectors for each
        %                       epoch in the return variable epoch, which
        %                       omits root.epoch for those shorter than
        %                       max_lag, and marges any overlapping epochs
        %                       in root.epoch so that we don't double count
        %
        % RETURNS
        %   cor                 either a vector (if one epoch, or average_epochs==1, or a column array where each column is 
        %                       the xcorr signal for epochs with spikes and
        %                       as long as max_lag
        %   lag                 same as above, but lags in seconds
        %   epoch               returns epochs used in analysis (in case
        %                       any of root.epoch were < max_lag, or had no
        %                       spikes
        %   std_error           if average_epochs=1, then the average error
        %                       for each lag timepoint is returned for the averaged cor vector
        %
        % andrew 3 april 2010
        
            import CMBHOME.* % has the spk_xcorr function
            import CMBHOME.Utils.*
            
            std_error = [];
            
            if nargin<2
                error('spk_xcorr: Not enough arguments');
            end
            
            if nargin<4
                t_bin=.001;
            end
            
            if nargin<3
                max_lag = 5; %seconds
            end
            
            if size(cells,2)~=2
                error('spk_xcorr: cells must be like [tetrode ind, cell ind]');
            end
            
            if size(cells,1)==1
                cells(2,:) = cells(1,:);
            end
            
            if ~exist('norm', 'var'), norm = 'count'; end

            if ~exist('average_epochs', 'var'), average_epochs = 0; end
            
            if ~exist('speed_dur', 'var'), speed_dur = 1; end
            
            if average_epochs, self = MergeEpochs(self); end % if we are averaging epochs, make sure that we aren't double counting
            
            epoch = self.epoch;
            
            tooshort = find( epoch(:,2)-epoch(:,1) < max_lag ); % remove any epochs which are shorter than max_lag
            
            self.epoch(tooshort, :) = []; epoch = self.epoch;

            if average_epochs, self = MergeEpochs(self); end % if we are averaging epochs, make sure that we aren't double counting
            
            if isempty(self.epoch)
                
                spk_xcorr = []; lags = []; epoch = []; std_error = [];
                
                warning('CMBH:error', 'No epochs meet this running speed requirement');
                
                return
                
            end
            
            x = self.spk_ts(cells(1,:));
            y = self.spk_ts(cells(2,:));
                    
            if iscell(x)
                    
                empty_inds = cellfun(@isempty, x); % find empty epochs and remove them

                empty_inds = empty_inds+cellfun(@isempty, y);

                empty_inds = logical(empty_inds);

                x(empty_inds) = []; % remove empty epochs

                y(empty_inds) = []; % remove empty epochs
                
                if isempty(x) % if no epochs with spikes
                    
                    lags = zeros(nbins, 1);
                    
                    spk_xcorr = zeros(nbins, 1);
                    
                    return
                    
                end                    

                epoch = epoch(~empty_inds, :); % remove epochs without spikes
                                                    
                nbins = length(-max_lag+t_bin/2:t_bin:max_lag-t_bin/2);
                
                spk_xcorr = zeros(nbins, length(x));
                
                lags = zeros(nbins, length(x));

                for i = 1 : length(x)
                    [spk_xcorr(:, i), lags(:, i)] = Spike.CrossCorr(x{i}, y{i}, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', norm);
                end
                
                if average_epochs % cross correlation across all epochs
                   
                    spk_xcorr = sum(spk_xcorr, 2);
                                       
                    lags = lags(:,1);

                    
                end
                
            else          
            
                [spk_xcorr, lags] = Spike.CrossCorr(x, y, 'lag', [-max_lag max_lag], 'binsize', t_bin, 'norm', norm);
                
            end

        end
        
        function vec_theta = spk_theta(self, cel)
        % spk_theta = root.spk_theta(cel);
        %
        % Returns spike theta phase of cel for every root.epoch. Vector if
        % one cell, and one root.epoch, cell array of size MxN where M is
        % the number of epochs and N is the number of cells
        %
        % cel -> [tetrode index, cell index]
        %
        % andrew 25 may 2010
        
        vec_theta = []; % initialize 
        
        ind_lfp = self.active_lfp;
                    
            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end
            
            if ~isempty(self.spike)
                
                if ~all([size(self.spike,1)>=cel(:,1), size(self.spike,2)>=cel(:,2)])
                    warning('CMBH:error', 'Indeces exceed number of tetrodes or cells');
                    return;
                end
                
            end
            
            if size(self.epoch,1)>1 || size(cel,1)>1
                vec_theta = cell(size(self.epoch,1) , size(cel,1));
            end
            
            for j = 1:size(cel,1)
                   
                [~, i] = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

                if iscell(i)
                    
                    i = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).i_lfp{ind_lfp}(c), i, 'Unif', 0);
                    
                    vec_theta(:,j) = cellfun(@(c) self.b_lfp(ind_lfp).theta_phase(c), i, 'UniformOutput', false);
                    
                elseif size(cel,1)>1
                    
                    i = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).i_lfp{ind_lfp}(c), i, 'Unif', 0);
                    
                    vec_theta{1,j} = self.b_lfp(ind_lfp).theta_phase(i);
                else % i is a cell array of indeces
                    
                    i = self.spike(cel(j,1), cel(j,2)).i_lfp{ind_lfp}(i);
                    
                    vec_theta = self.b_lfp(ind_lfp).theta_phase(i);
                end

            end
            
        end
        
        function vec_x = spk_x(self, cel, legacy)
        % spk_x = root.spk_x(cel);
        %
        % Returns spike x position of cel for every root.epoch. Vector if
        % one cell, and one root.epoch, cell array of size MxN where M is
        % the number of epochs and N is the number of cells
        %
        % cel -> [tetrode index, cell index]
        %
        % andrew 25 may 2010
        
        vec_x = []; % initialize 
        
            if nargin==3
                cel = [cel, legacy];
                help CMBHOME.Session.spk_x
                warning('CMBH:error', 'Your command is outdated to match toolbox convention');
            end
            
            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end
            
            if ~isempty(self.spike)
                
                if ~all([size(self.spike,1)>=cel(:,1), size(self.spike,2)>=cel(:,2)])
                    warning('CMBH:error', 'Indeces exceed number of tetrodes or cells');
                    return;
                end
                
            end
            
            if size(self.epoch,1)>1 || size(cel,1)>1
                vec_x = cell(size(self.epoch,1) , size(cel,1));
            end
            
            for j = 1:size(cel,1)
            
                i = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

                if iscell(i)
                    vec_x(:,j) = cellfun(@(c) self.b_x(c), i, 'UniformOutput', false);
                elseif size(cel,1)>1
                    vec_x{1,j} = self.b_x(i);
                else % i is a cell array of indeces
                    vec_x = self.b_x(i);
                end

            end
            
        end
        
        function vec_y = spk_y(self, cel, legacy)
        % spk_y = root.spk_y(cel);
        %
        % Returns spike y position of cel for every root.epoch. Vector if
        % one cell, and one root.epoch, cell array of size MxN where M is
        % the number of epochs and N is the number of cells
        %
        % cel -> [tetrode index, cell index]
        %
        % andrew 25 may 2010
        
        vec_y = []; % initialize
        
            if nargin==3
                cel = [cel, legacy];
                help CMBHOME.Session.spk_y
                warning('CMBH:error', 'Your command is outdated to match toolbox convention');
            end
            
            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end
            
            if ~isempty(self.spike)
                
                if ~all([size(self.spike,1)>=cel(:,1), size(self.spike,2)>=cel(:,2)])
                    warning('CMBH:error', 'Indeces exceed number of tetrodes or cells');
                    return;
                end
                
            end
            
            if size(self.epoch,1)>1 || size(cel,1)>1
                vec_y = cell(size(self.epoch,1) , size(cel,1));
            end
            
            for j = 1:size(cel,1)
            
                i = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

                if iscell(i)
                    vec_y(:,j) = cellfun(@(c) self.b_y(c), i, 'UniformOutput', false);
                elseif size(cel,1)>1
                    vec_y{1,j} = self.b_y(i);
                else % i is a cell array of indeces
                    vec_y = self.b_y(i);                
                end
                
            end
            
        end
        
        function spk_i = spk_i(self, cel)
           % spk_ts = root.spk_ts(cel);
        %
        % Returns spike timestamps of cel for every root.epoch. Vector if
        % one cell, and one root.epoch, cell array of size MxN where M is
        % the number of epochs and N is the number of cells
        %
        % cel -> [tetrode index, cell index]
        %
        % andrew 25 may 2010
        
            spk_i = []; % initialize
           
            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end
            
            if ~isempty(self.spike)
                
                if ~all([size(self.spike,1)>=cel(:,1), size(self.spike,2)>=cel(:,2)])
                    warning('CMBH:error', 'Indeces exceed number of tetrodes or cells');
                    return;
                end
                
            end
            
            if size(self.epoch,1)>1
                spk_i = cell(size(self.epoch,1) , size(cel,1));
            end
            
            for j = 1:size(cel,1)
                
                i = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

                if iscell(i)
                    spk_i(:,j) = i;
                elseif size(cel,1)>1
                    spk_i{1,j} = i;
                else % i is a vector and there is only one cell
                    spk_i = i;
                end

            end   
            
        end
        
        function vec_t = spk_ts(self, cel, legacy)
        % spk_ts = root.spk_ts(cel);
        %
        % Returns spike timestamps of cel for every root.epoch. Vector if
        % one cell, and one root.epoch, cell array of size MxN where M is
        % the number of epochs and N is the number of cells.
        %
        % Note that some spikes may exist in root.epoch, but due to
        % downsampling to root.fs_video, now exist outside of the epoch
        % period. Those spikes are not included, so that all of the spk_*
        % functions return vectors of the same length.
        %
        % cel -> [tetrode index, cell index]
        %
        % andrew 25 may 2010
        
        vec_t = []; % initialize
        
            if nargin==3
                cel = [cel, legacy];
                help CMBHOME.Session.spk_ts
                warning('CMBH:error', 'Your command is outdated to match toolbox convention');
            end
            
            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end
            
            if ~isempty(self.spike)
                
                if ~all([size(self.spike,1)>=cel(:,1), size(self.spike,2)>=cel(:,2)])
                    warning('CMBH:error', 'Indeces exceed number of tetrodes or cells');
                    return;
                end
                
            end
            
            if size(self.epoch,1)>1
                vec_t = cell(size(self.epoch,1) , size(cel,1));
            end
            
            for j = 1:size(cel,1)

                [~, i] = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);
                %i = IsolateEpoch(self.spike(cel(j,1), cel(j,2)).ts, 1:length(self.spike(cel(j,1), cel(j,2)).ts), self.epoch);

                if iscell(i)
                    vec_t(:,j) = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).ts(c), i, 'UniformOutput', false);
                elseif size(cel,1)>1
                    vec_t{1,j} = self.spike(cel(j,1), cel(j,2)).ts(i);
                else % i is a vector and there is only one cell
                    vec_t = self.spike(cel(j,1), cel(j,2)).ts(i);
                end

            end
            
        end
        
        function vec_bin = spk_bin(self, cel, dt)
        % spike_train_binned = root.spk_bin(cel, dt);
        %
        % Returns binned spike observations for all root.epoch-s at
        % binwidth dt (default 2ms). If multiple epochs, vec_bin is a cell
        % array
        %
        % cel -> 2 element vector like [tetrode index, cell index]
        % dt -> positive time binwidth
        %
        % andrew 25 may 2010
        
            import CMBHOME.*

            if ~exist('cel', 'var'), 
                help CMBHOME.Session.spk_bin
                error('You must indicate a cell.'); 
            end

            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end

            if ~exist('dt', 'var'), dt = .002; end

            epochs = self.epoch;
            
            vec_bin = cell(size(epochs,1), 1);

            for i = 1:size(epochs,1)

                self.epoch = epochs(i,:); % set epoch to what we're interested in

                vec_bin{i} = Spike.TS2Bins(self.spk_ts(cel), self.epoch, dt);

            end

            if length(vec_bin)==1, vec_bin = vec_bin{1}; end % convert to array if only one epoch
        
            
        end
        
        function vec_vel = spk_vel(self, cel)
        % spk_vel = root.spk_vel(cel);
        %
        % Returns spike running speed of cel for every root.epoch. Vector if
        % one cell, and one epoch, cell array of size MxN where M is
        % the number of epochs and N is the number of cells
        %
        % cel -> [tetrode index, cell index]
        %
        % andrew 25 may 2010
        
            vec_vel = []; %initialize

            if ~isempty(self.b_vel), vel = self.b_vel; 
            else vel = cat(1, 0, sqrt(diff(self.b_x).^2+diff(self.b_y).^2))/self.fs_video^-1;
            end

            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end

            if ~isempty(self.spike)

                if ~all([size(self.spike,1)>=cel(:,1), size(self.spike,2)>=cel(:,2)])
                    warning('CMBH:error', 'Indeces exceed number of tetrodes or cells');
                    return;
                end

            end

            if size(self.epoch,1)>1
                vec_vel = cell(size(self.epoch, 1) , size(cel,1));
            end
            
            for j = 1:size(cel,1)
            
                i = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

                if iscell(i)
                    vec_vel(:,j) = cellfun(@(c) vel(c), i, 'UniformOutput', false);
                elseif size(cel,1)>1
                    vec_vel{1,j} = vel(i);
                else % i is a cell array of indeces               
                    vec_vel = vel(i);                
                end
                
            end
            
        end
        
        function vec_headdir = spk_headdir(self, cel, legacy)
        % spk_headdir = root.spk_y(cel);
        % OLD SYNTAX: root.spk_headdir(tetrode_i, cell_i)
        %
        % Returns spike y position of cel for every root.epoch. Vector if
        % one cell, and one root.epoch, cell array of size MxN where M is
        % the number of epochs and N is the number of cells
        %
        % Arguments:
        % cel -> Nx2 vector indicating [tetrode, cell].
        %
        % requires root.epoch -> Mx2 array of epoch times.
        %
        % Returns:
        %
        % headdir -> vector if M and N are 1, cell array of size MxN
        % otherwise (epochs x cells)
        %
        % andrew 25 may 2010
        
            vec_headdir = []; % initialize
            
            if isempty(self.b_headdir), error('No head direction data'); end
        
            if nargin==3
                cel = [cel, legacy];
                help CMBHOME.Session.spk_headdir
                warning('CMBH:error', 'Your command is outdated to match toolbox convention');
            end
            
            if numel(cel)<2, error('Cell must be implicated by tetrode index and cell index.'); end
            
            if ~isempty(self.spike)
                
                if ~all([size(self.spike,1)>=cel(:,1), size(self.spike,2)>=cel(:,2)])
                    warning('CMBH:error', 'Indices exceed number of tetrodes or cells');
                    return;
                end
                
            end
            
            if size(self.epoch,1)>1 || size(cel,1)>1
                vec_headdir = cell(size(self.epoch,1) , size(cel,1));
            end
            
            for j = 1:size(cel,1)

                i = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

                if iscell(i)
                    vec_headdir(:,j) = cellfun(@(c) self.b_headdir(c), i, 'UniformOutput', false);
                elseif size(cel,1)>1
                    vec_headdir{1,j} = self.b_headdir(i);
                else % i is a cell array of indeces
                    vec_headdir = self.b_headdir(i);
                end

            end
                
        end
        
        function self = AppendKalmanVel(self)
            
            if ~isempty(self.b_vel)
                
                warning('CMBH:error', 'User defined speed vector already created. Clear root.b_vel first.');
                
                return
                
            else
            
                [t,x,y,vx,vy] = self.KalmanVel(self.b_x,self.b_y,self.b_ts, 2); % quadratic

                [trash, inds] = ismember(t, self.b_ts); % all indices in self.b_ts for which we have a kalman estimate

                b_vel = zeros(size(self.b_ts));

                b_vel(inds) = sqrt(vx.^2 + vy.^2); % pixels/sec
                
                self.b_vel = b_vel;
                
                warning('CMBH:notify', '%s', 'Added Kalman Velocity to Session Object.');
            
            end
        end
        
        function ts = Label2Time(self, label)
        % ts = root.Label2Time(str)
        %
        % Returns ts, a vector of timestamp(s) from root.event that line up
        % with event label 'str'
        
        ts = [];
        
        if iscell(label)
            
            ts = zeros(size(label, 1), size(label,2));
            
            for i = 1:size(label, 1)
                for j = 1:size(label, 2)
                    tmp = cell2mat(self.event(ismember(self.event(:,1), label{i,j}),2));
                    
                    if numel(tmp)==1, ts(i,j) = tmp;
                    elseif isempty(tmp), ts(i, j) = NaN;
                    else error('Multiple identical labels exist, cannot use cell array functionality');
                    end
                    
                end
            end
        else
        
            ts = cell2mat(self.event(ismember(self.event(:,1), label),2));
        
        end

        end
        
        function self = AddLFP(self, fname)
        % (1) root = root.AddLFP
        % (2) root = root.AddLFP(fname)
        %
        % Appends root.path_lfp with valid LFP file name via user input
        % prompts or filename, and adds it to root.path_lfp.
        %
        % (1) Prompts the user to select a file using the file browser
        % (2) Checks if string 'fname' is a valid LFP file (.ncs, .mat,
        % .plx files) and adds it to root.path_lfp.
        %
        % andrew 3 april 2010
                        
            if ~exist('fname', 'var')
                if isempty(self.path_lfp) 
                    i = 1;
                else
                    i = length(self.path_lfp)+1;
                end

                tmppath = self.name;

                tmppath = tmppath(1:find(tmppath==filesep, 1, 'last')-1); % directory of datafile

                [load_file,base_path] = uigetfile('*.eeg*;*.ncs;*.mat;*.plx','Select the LFP (EEG) file relevent to this project. Cancel to exit',tmppath, 'MultiSelect', 'on');

                if iscellstr(load_file)
                    
                    for j = 1:length(load_file)
                        self.path_lfp{i} = fullfile(base_path, load_file{j});
                        i = i+1;
                    end

                    self = self.AddLFP; % recursively call this function to load files
                    
                    warning('CMBH:notify', '%s', 'Added'), warning('CMBH:notify', '%s', load_file(:));
                    
                elseif load_file
                                       
                    self.path_lfp{i} = fullfile(base_path, load_file);

                    self = self.AddLFP; % recursively call this function to load files
                    
                    warning('CMBH:notify', '%s', 'Added'), warning('CMBH:notify', '%s', load_file(:));
                    
                end
            else
                
                if ~iscell(fname), fname = {fname}; end % convert to cell if char
                
                for i = 1:length(fname)
                    
                    if exist(fname{i}, 'file')
                        self.path_lfp{end+1} = fname{i};
                    else
                        warning('CMBH:error', [fname{i} ' does not exist. Not added']);
                        fname{i} = '';
                    end
                    
                end
                
                warning('CMBH:notify', '%s', ['Added ' int2str(length(fname)) ' lfp filenames.']);
                
            end
                                       
        end  
        
        function self = ClearLFP(self)
        % root = root.ClearLFP;
        %
        % Clears the LFP signals from root.b_lfp. This may be useful if the
        % LFP signals are very large, or highly sampled.
        %
        % andrew 3 april 2010
        
            self.b_lfp = [];
            
            self.active_lfp = [];
                         
        end 
        
        function Status(self)
        % root.Status
        % prints to screen relevant information about this Session class
        % object
        %
        % andrew 3 april 2010

            switch self.cell_thresh(2)
                case 1
                    str_cell = ' Hz continuously throughout the entire session';
                case 0
                    str_cell = ' Hz on average during the session';
            end

            disp(['This CMB Session is ' self.name '.' char(10) 'It has ' int2str(size(self.cells,1)) ...
                ' cells which fire above ' num2str(self.cell_thresh(1)) ...
                str_cell '.' char(10) char(10) 'There are ' int2str(size(self.event,1)) ' event flags' ...
                '.' char(10) 'The recording duration is ' num2str(self.b_ts(end)-self.b_ts(1)) ' seconds.']);
            
            disp([char(10) 'It was created on ' datestr(self.date_created, 1)]);
            
            disp([char(10) 'The tracking sampling frequency is ' num2str(self.fs_video) ' Hz.']);
            
            if isempty(self.b_vel), disp('root.vel returns the simple derivative of position. Edit root.b_vel to specify a more accurate vector of running speed,');
            else, disp('User defined speed vector'); end
            
            switch self.raw_headdir
                case 1
                    str_headdir = 'The head direction signal has not been corrected or smoothed';
                case 0
                    str_headdir = 'The head direction signal has been corrected and smoothed';
            end
            
            switch self.raw_pos
                case 1
                    str_pos = 'The position data has not been corrected or smoothed.';
                case 0
                    str_pos = 'The position data has been corrected and smoothed.';
            end
            
            switch isempty(self.spatial_scale)
                case 1
                    str_spat_scale = [char(10) 'You have not indicated resolution (cm/pixel). ' ...
                        'Do so by assigning object.spatial_scale a value.'];
                case 0
                    str_spat_scale = [char(10) 'The resolution is set to ' num2str(self.spatial_scale) ...
                        ' cm/pixel.'];
            end
            
            disp(str_spat_scale);            
            disp(str_headdir);
            disp(str_pos);        
            
        end
        
        function self = SetEpoch(self)
        % root = root.SetEpoch;
        %
        % GUI for setting epoch start and stop times for plotting, etc
           
            width = 425;
            height = 600;

            %  Create and then hide the GUI as it is being constructed.
            f = figure('Visible','off','Position', [50, -150, width, height],'Color', 'w');

            sizes = [15, 180, 260];
            widths = [150, 75, 150];
            lineheight = 30;
            
            check_ind = 1;

            %  Construct the components.
            
            %labels
            hlabel1 = uicontrol('Style','text','String','Anchoring Event',...
                  'Position',[sizes(1),height-50,widths(1),lineheight], 'BackgroundColor', 'w',...
                  'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
            hlabel2 = uicontrol('Style','text','String','Lag (+/- seconds)',...
              'Position',[sizes(2),height-50,widths(2),lineheight], 'BackgroundColor', 'w',...
              'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
            hlabel3 = uicontrol('Style','text','String','Secondary Event',...
              'Position',[sizes(3),height-50,widths(3),lineheight], 'BackgroundColor', 'w',...
              'FontSize', 8, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
              
            str_help = ['If you specify start and stop above, we will ignore the selected labels.' ...
                ' Start and stop must lie between 0 and 1, and indicate percent session time.'];
          
            helplabel = uicontrol('Style','text','String',str_help,...
              'Position',[5,65,width-10,60], 'BackgroundColor', 'w',...
              'FontSize', 10, 'HorizontalAlignment', 'left', 'HandleVisibility', 'off');
          
            h_objarr(check_ind, 1) = uicontrol('Style','popupmenu',... % one popup
                  'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)),char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                  'Position',[sizes(1),height-50-check_ind*lineheight,widths(1),lineheight],...
                  'HandleVisibility', 'off');
              
            h_objarr(check_ind, 2) = uicontrol('Style','edit',... % one textbox
                  'String','',...
                  'Position',[sizes(2),height-50-check_ind*lineheight+8,widths(2),22],...
                  'HandleVisibility', 'off');

            h_objarr(check_ind, 3) = uicontrol('Style','popupmenu',... % one popup
                  'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)),char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                  'Position',[sizes(3)+15,height-50-check_ind*lineheight,widths(3), lineheight],...
                  'HandleVisibility', 'off');
         
            h_objarr(check_ind, 4) = uicontrol('Style','edit',... % one textbox
                  'String','Start',...
                  'Position',[sizes(1),height-50-check_ind*lineheight*1.5,widths(2),22],...
                  'HandleVisibility', 'off');  
              
            h_objarr(check_ind, 5) = uicontrol('Style','edit',... % one textbox
                  'String','Stop',...
                  'Position',[sizes(2)-50,height-50-check_ind*lineheight*1.5,widths(2),22],...
                  'HandleVisibility', 'off');

              
            % buttons for add epoch or finish
            
            h_add = uicontrol('Style','pushbutton','String','Add',...
                  'Position',[width-100, height-50-(check_ind+1)*lineheight*1.5, 50, 25],...
                  'Callback',{@add_Callback}, 'HandleVisibility', 'off'); 
            h_finish = uicontrol('Style','pushbutton','String','Done',...
                  'Position',[width-50, height-50-(check_ind+1)*lineheight*1.5, 50, 25],...
                  'Callback',{@finish_Callback}, 'HandleVisibility', 'off');

            % Initialize the GUI.
            % Change units to normalized so components resize 
            % automatically.

            % set([f; h_objarr; h_add; h_finish; hlabel1; hlabel2; hlabel3],'Units','Pixels');

            %ha = axes('Position', [.1 .1 .6 .6]);

            % Assign the GUI a name to appear in the window title.
            set(f,'Name','SetEpoch')
            % Move the GUI to the center of the screen.
            movegui(f,'center')
            % Make the GUI visible.
            set(f,'Visible','on');
            
            uiwait(f); % wait to return until figure is closed
            
            function popup_lag_Callback(source, eventdata) 
          
                % Determine the selected data set.
                % clear old plots, and update header

                str = get(source, 'String');
                val = get(source,'Value');
                
            end
            
            function add_Callback(source, eventdata) % add another line for epoch
                
                check_ind = check_ind+1;
                
                h_objarr(check_ind, 1) = uicontrol('Style','popupmenu',... % one popup
                      'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)), char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                      'Position',[sizes(1),height-15-check_ind*lineheight*2,widths(1), lineheight],...
                      'HandleVisibility', 'off');
              
                h_objarr(check_ind, 2) = uicontrol('Style','edit',... % one textbox
                      'String','',...
                      'Position',[sizes(2),height-15-check_ind*lineheight*2+8,widths(2),22],...
                      'HandleVisibility', 'off');

                h_objarr(check_ind, 3) = uicontrol('Style','popupmenu',... % one popup
                      'String',cellstr(char(cat(2, num2str(shiftdim(1:size(self.event,1),1)),char(self.event{:,1}), num2str(vertcat(self.event{:,2}))))),...
                      'Position',[sizes(3),height-15-check_ind*lineheight*2,widths(3), lineheight],...
                      'HandleVisibility', 'off');
                  
                h_objarr(check_ind, 4) = uicontrol('Style','edit',... % one textbox
                  'String','Start',...
                  'Position',[sizes(1),height-15-check_ind*lineheight*2-lineheight+8,widths(2),22],...
                  'HandleVisibility', 'off');  
              
                h_objarr(check_ind, 5) = uicontrol('Style','edit',... % one textbox
                  'String','Stop',...
                  'Position',[sizes(2)-50,height-15-check_ind*lineheight*2-lineheight+8,widths(2),22],...
                  'HandleVisibility', 'off');
                  
                set(h_add, 'Position', [width-100, height-15-(check_ind)*lineheight*2 - lineheight+8, 50, 25]);
                set(h_finish, 'Position', [width-50, height-15-(check_ind)*lineheight*2 - lineheight+8, 50, 25]);   
                
            end
            
            function finish_Callback(source, eventdata) % set self.epoch and close window
                
                epoch_ind = 1;
                
                dur = self.b_ts(end)-self.b_ts(1); % duration of recording

                for i = 1:size(h_objarr,1)
                    
                    user_start = get(h_objarr(i,4), 'String');
                    user_stop = get(h_objarr(i,5), 'String'); 
                    
                    if ~strcmp(user_start, 'Start') && ~strcmp(user_stop, 'Stop') %start and stop are declared and valid, 
                        
                        user_start = str2num(user_start);
                        user_stop = str2num(user_stop);
                        
                        if user_start<user_stop && user_start<=1 && user_start>=0 && user_stop<=1 && user_stop>=0
                            epoch(epoch_ind,:) = [dur*user_start+self.b_ts(1), dur*user_stop+self.b_ts(1)];
    
                            epoch_ind = epoch_ind+1;
                        else
                            disp('Start and Stop must lie between 0 and 1. Defaulting to 0 and 1 (start and end of session)');
                        end
                        
                    else
                    
                           
                        t1ind = get(h_objarr(i,1), 'Value');
                        t1 = self.event{t1ind,2};

                        if ~isempty(get(h_objarr(i,2), 'String'))
                            t2 = t1+str2num(get(h_objarr(i,2), 'String'));
                            if t2<0
                                disp('One of the timestamps was less than zero. Setting to zero');
                                t2=0;
                            end
                        end

                        if ~exist('t2', 'var')                        
                            t2ind = get(h_objarr(i,3), 'Value');
                            t2 = self.event{t2ind,2};
                            clear t2ind
                        end

                        if t1~=t2
                            epoch(epoch_ind,1) = min([t1, t2]);
                            epoch(epoch_ind,2) = max([t1, t2]);
                            epoch_ind = epoch_ind+1;
                        end
                        clear t1 t2 ind 
                        
                    end
                end
                
                self.epoch = epoch;
                
                close(f)
                
            end
            
        end
        
        function self = LoadLFP(self, ind)
        % (1) root = root.LoadLFP
        % (2) root = root.LoadLFP(indices)
        %
        % Loads LFP data file into root.b_lfp property
        %
        % (1) Loads root.path_lfp(root.active_lfp)
        % (2) Loads all root.path_lfp(indices) to current object
        %
        % So far, .mat files with lfp, lfp_ts fields can be imported. .ncs 
        % files can be imported. other methods
        % can be addded
        %
        % andrew 3 april 2010
        
            if ~exist('ind', 'var')
                
                if isempty(self.active_lfp), error('root.LoadLFP: Please specify root.active_lfp'); 
                else ind = self.active_lfp; end
                
            elseif isempty(ind)
                
                warning('CMBH:error', 'root.LoadLFP: Indeces in root = root.LoadLFP(indeces) was empty, so nothing was done. Do not include argument indeces to use root.active_lfp.');
                
                return
            
            end
            
            
        
            if ~iscell(self.path_lfp), error('root.LoadLFP: Please specify LFP files via root.AddLFP'); end
            
            for i = 1:length(ind)
            
                f_name = self.path_lfp{ind(i)};

                if ~strcmp(class(self.b_lfp), 'CMBHOME.LFP'), self.b_lfp = CMBHOME.LFP; end

                if exist(f_name, 'file')       

                    if strfind(f_name(find(f_name=='.', 1, 'last'):end), '.eeg') % ehrens misbehaved data names
                        
                        % Read EEG file
                        %
                        % [EEG,Fs] = geteeg(datafile);
                        %
                        % You can create a time vector for plotting your EEG by taking 
                        % t = (0:length(EEG)-1)'/Fs
                        %
                        %
                        % Sturla Molden <sturla@molden.net>
                        % Centre for the Biology of Memory
                        % Norwegian University of Science and Technology
                        % http://www.cbm.ntnu.no
                        % 
                        % Copyright (C) 2003  Centre for the Biology of Memory, NTNU
                        % All Rights Reserved
                        %
                        % This M-file is released under Q Public License v 1.0,
                        % a copy of which should accompany this file.

                        fid = fopen(f_name,'r');
                        if (fid == -1)
                           error(sprintf('Could not open file %s',f_name));
                        end
                        for ii = 1:8
                           textstring = fgetl(fid);
                        end
                        fs = sscanf(textstring,'%*s %f');
                        for ii = 1:3
                           textstring = fgetl(fid);
                        end
                        nosamples = sscanf(textstring,'%*s %u');
                        fseek(fid,10,0);
                        signal = fread(fid,nosamples,'int8');
                        ts = (0:length(signal)-1)'/fs;
                        fclose(fid);

                        channel_name = f_name;

                        % theta_phase = DATA(:,7);
                        
                    else % some other, more reasonable 3 character extension
                    
                        switch lower(f_name(find(f_name=='.', 1, 'last'):end))
                            
                            case '.mat'

                                load(f_name)

                                channel_name = ind(i);

                                if exist('lfp', 'var')
                                    signal = lfp;
                                else
                                    warning('CMBH:error', 'Vector with signal not found.')
                                end

                                if exist('lfp_ts', 'var')
                                    ts = lfp_ts;
                                    fs = mean(diff(lfp_ts))^-1;
                                else
                                    warning('CMBH:error', 'Vector with signal timestamps not found. Cannot assign root.fs')
                                end


                            case '.plx'
                                % not written yet
                                warning('CMBH:error', 'Not coded up yet');

                            case '.ncs'
                                % wow, this neuralynx literature is a mess. ok:
                                % 'n_samples' below is a vector with an integer indicating
                                % how many samples exist at a certain index that
                                % you can pass at the fifth parameter to NlxMatCSC
                                import CMBHOME.NL2Mat.*

                                import_segment_length = 2^14; % each segment has 2^12 samples

                                [fs, n_samples] = Nlx2MatCSC(f_name, [0 0 1 1 0], 0, 1); % check the sampling frequency, and get the number of samples
                                %     1. Timestamps   
                                %     2. Channel Numbers
                                %     3. Sample Frequency
                                %     4. Number of Valid Samples
                                %     5. Samples

                                if length(unique(fs))~=1, warning('CMBH:notify', '%s', ['Sampling rate in ' fname ' is not consistent. Check timestamps.']); end  

                                decimation_factor = round(fs(1)/400); % we want to downsample to around 400 Hz

                                decimation_factor = 2^(nextpow2(decimation_factor)-1); % simplifies our downsampling between time and signal

                                fs_new = fs(1)/decimation_factor;

                                if length(n_samples)*512 ~= sum(n_samples), warning('CMBH:notify', '%s', 'Tell Andrew Neuralynx is doing something funny.'); end

        %                         signal = zeros(length(n_samples) * 512, 1);
        %                         tmp_ts = zeros(length(n_samples), 1);

                                %signal = zeros(ceil((512 * import_segment_length)/decimation_factor) * floor(n_samples/import_segment_length)... % build signal vector
                                   % + ceil(((512 * import_segment_length) - floor(n_samples/import_segment_length) * 512)/decimation_factor), 1);

                                warning('CMBH:notify', '%s', ['Downsampled from ' num2str(fs(1)) 'Hz to ' num2str(fs_new) 'Hz.']); 

                                [tmp_ts, samples] = Nlx2MatCSC(f_name, [1 0 0 0 1], 0, 1);

        %                         i = 0;
        %                         while i < length(n_samples) % extract signal in chunks
        %                             
        %                             inc = min([length(n_samples) - i, import_segment_length]);
        %                                 
        %                             [ttmp_ts, tmp_samples] = Nlx2MatCSC(f_name, [1 0 0 0 1], 0, 2, [i i+inc-1]);
        %                             
        %                             tmp_ts(i+1:i+inc) = ttmp_ts;
        %                             
        %                             signal(i*512+1 : (i+inc)*512) = tmp_samples(:);
        %                             
        %                             i = i+inc;
        %                             
        %                         end

                                signal = downsample(samples(:), decimation_factor);

                                l_seg = floor(length(signal)/length(tmp_ts));

                                ts = zeros( length(tmp_ts) * l_seg , 1 );

                                tmp_ts(end+1) = tmp_ts(end)+mean(diff(tmp_ts));

                                for ii = 1:length(tmp_ts)-1

                                    dt = (tmp_ts(ii+1)-tmp_ts(ii))/l_seg * (0:l_seg);
                                    dt = dt(1:end-1);
                                    ts( (ii-1) * l_seg + 1 : ii * l_seg) = tmp_ts(ii) + dt;

                                end     

        %                         indices = 0:import_segment_length:length(n_samples)-1; % build indices
        %                         
        %                         samples = zeros(floor(sum(n_samples)/decimation_factor),1);
        %                         ts = zeros(size(samples));
        %                         
        %                         if indices(end)~=length(n_samples)-1, indices(end+1) = length(n_samples)-1; end
        %                         
        %                         for i = 1:length(indices)-1 % parse through all indeces
        %                              [tmp_ts, tmp_samples] = Nlx2MatCSC(f_name, [1 0 0 0 1], 0, 2, [indices(i) indices(i+1)-1]);
        % 
        %                              l_seg = floor(import_segment_length / decimation_factor);
        %                              
        %                              tmp_samples = downsample(tmp_samples(:), decimation_factor);
        %                                                           
        %                              samples(sum(indices(1:i))+1:sum(indices(1:i))+l_seg) = tmp_samples;
        %                              ts(sum(indices(1:i))+1:sum(indices(1:i))+l_seg) = tmp_ts(1)+(0:l_seg-1)*fs_new^-1;                    
        %                         end

                                ts = ts*10^-6; % convert from microseconds to seconds

        %                         for i = 1:length(ts)-1
        %                             lfp_ts(1+(i-1)*fold:i*fold) = ts(i) : (ts(i+1)-ts(i))/(fold-1) : ts(i+1); 
        %                         end

                               % reshape the lfp signal
                                ts = ts(:);
                                channel_name = self.path_lfp{ind(i)};
                                fs = fs_new;

                        end
                    end
                    
                    self.b_lfp(ind(i)) = CMBHOME.LFP(signal, ts, fs, channel_name);
                    

                else
                    warning('CMBH:error', ['It appears file ' self.path_lfp{ind} ' in root.path_lfp does not exist, and was skipped. Run root.AddLFP to find LFP file(s) again.']);
                end
            
            end
            
        end
        
        function Save(root, params) 
        % (1) root.Save
        % (2) root.Save(params)
        %
        % Saves Session class object to root.name_formatted
        %
        % (1) Saves root to root.name, asks if overwrite is necessary. If
        % LFP exists, asks if it should be cleared.
        %
        % (2) Saves root to root.name, but will not prompt user as per
        % params:
        %
        % params = [ <clear_lfp> , <noprompt> ] for batch processing ex.
        % root.Save([1 0]) will clear the LFP, but not overwrite root.name,
        % if it exists
        %
        % If noprompt == 1, then root.Save will attempt to save to
        % root.name. If Directory does not exist, it will create it. If
        % file exists, it will overwrite.
        %
        % If you want to save to a new filename, change root.name before
        % calling root.save
        
            if ~exist('params', 'var')
                params = [NaN NaN];
            elseif numel(params)==1
                params(2) = NaN;
            end
        
        
            if ~isempty(root.b_lfp)
                
                if params(1)==1
                    root = root.ClearLFP;
                    warning('CMBH:notify', '%s', 'Cleared LFP from object. If you want to save it in the future, use root.Save(0)');
                elseif params(2)~=1
                    user_entry = input('Would you like to save the LFP data with the object? This could make for a large filesize. [y/n]: ', 's');
                
                    while( ~strcmpi(user_entry, 'y') && ~strcmpi(user_entry, 'n') )
                        user_entry = input('Would you like to save the LFP data with the object? This could make for a large filesize. [y/n]: ', 's');
                    end
                    
                    if strcmpi(user_entry, 'n')
                        root = root.ClearLFP;
                        warning('CMBH:notify', '%s', 'Cleared LFP from object.');
                    end
                end
                
            end
                        
           
            if exist(root.name_formatted, 'file') && (isnan(params(2)) || params(2)==0) % if exists and user didnt specify what to do
                reply = input(['You are about to overwrite ' root.name_formatted '. Are you sure? [y/n]: '], 's');
                if strcmpi(reply, 'y')            
                    save(root.name_formatted, 'root', '-v7.3');
                    warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
                else
                    warning('CMBH:notify', '%s', 'Save aborted.');
                end               
            elseif exist(root.name_formatted, 'file') && params(2)==1 % if exists and user specified overwrite
                save(root.name_formatted, 'root', '-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
            elseif exist(root.name_formatted(1:find(root.name_formatted==filesep, 1, 'last')-1), 'dir') % if directory exists, but file doesnt, save it
                save(root.name_formatted, 'root', '-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
            elseif params(2)~=1 % if directory doesnt exist and noprompt is false, ask user
                disp('Folder in root.name not found. Please select where you would like to save the file.');
                [file, path] = uiputfile('*.mat', 'Save Object As', pwd);
                root.name = fullfile(path, file);
        
                save(root.name_formatted, 'root','-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
            else % directory doesnt exist, but noprompt is true, so mkdir and save file
                
                mkdir(root.name_formatted(1:find(root.name_formatted==filesep, 1, 'last')-1));
        
                save(root.name_formatted, 'root', '-v7.3');
                warning('CMBH:notify', '%s', ['Saved ' root.name_formatted]);
                
            end
            
        end
        
        function Session_array = SplitSession(self, varargin)
        % (1) Session_array = root.SplitSession;
        % (2) root.SplitSession(params);
        %
        % This method will split the present Session instance into more
        % than one Session Objects.
        %
        % (1) Splits the session object as per root.epochs and optionally
        % returns 'Session_array' of session objects of length(size(root.epoch,1))
        %
        % (2) Takes optional parameters
        % 'filenames' -  cellarray of strings, filenames to make each root.name
        % 'epoch' - Nx2 array of times to split. N = length off 'filenames'
        % 'save' - 1 or 0 (default 0). Save each session object
        %
        % andrew 24 may 2010
               
        import CMBHOME.*
        
        success = 0; % initialize variables
        Session_array = Session;
        
        p = inputParser;

        p.addRequired('self');
        p.addParamValue('filenames', [], @(x) isccellstr(x)); 
        p.addParamValue('epoch',  [], @(x) size(x,2)==2);
        p.addParamValue('save', 0, @ischar);

        p.parse(self, varargin{:});

        filenames = p.Results.filenames;
        epoch = p.Results.epoch;
        save = p.Results.save;
        
        if isempty(epoch)
            epoch = self.epoch;
        end
        
        if size(epoch,1)<2 % if there arent any epochs passed
            help CMBHOME.Session.SplitSession
            return;
        end
        
        if isempty(filenames)
            
            append = 1:size(epoch,1);
            
            filenames = cat(2, repmat(strrep(self.name, '.mat', ''), size(epoch,1), 1), repmat('_Split', size(epoch,1), 1), int2str(append'), repmat('.mat', size(epoch,1), 1));
            
            filenames = mat2cell(filenames, ones(size(filenames,1),1), size(filenames,2));

        end
        
        if length(filenames)~=size(epoch,1)
            help CMBHOME.Session.SplitSession
            error('length of filenames does not match number of epochs');
        end
        
        for i = 1:length(filenames) % parse through all Session Objects
            
            self.epoch = epoch(i,:);
            
            self.cell_thresh = [ 0 0 ];
            
            event = self.event( [self.event{:,2}]>= epoch(i,1) & [self.event{:,2}]<=epoch(i,2), : );
        
            Session_array(i) = Session('name', filenames{i},...
                                    'b_ts', self.ts,'b_x',self.x, 'b_y', self.y, 'b_headdir',...
                                    self.headdir, 'event', event, 'raw_pos', self.raw_pos, ...
                                    'raw_headdir', self.raw_headdir, 'date_created', now, ...
                                    'epoch', epoch(i,:), 'fs_video', self.fs_video,...
                                    'path_lfp', self.path_lfp, 'path_raw_data', self.path_raw_data); % put it all together

            Session_array(i).user_def = self.user_def; % copy over all user_def fields                    
                                
            Session_array(i).user_def.SplitHistory = ['This Session object was split from ' self.name ' using SplitSession.'];                   
                                
            for j = 1:size(self.cells,1) % import cells
                
                spike = self.spike(self.cells(j,1), self.cells(j,2));
                
                spike.ts = spike.ts(spike.ts>=epoch(i,1) & spike.ts<=epoch(i,2));
                
                spike.i = Utils.SpkInds(self.ts, spike.ts);
               
                Session_array(i).spike(self.cells(j,1), self.cells(j,2)) = spike;
                
            end
            
            if save, Session_array(i).Save; end % save, if specified                                
            
        end           
            
    end
    end
    
    properties (Hidden) % old or outdated properties, around for backwards compatibility
    
            b_kalmanvel
            
    end
end

function [i_ts, i_i] = IsolateEpoch(ts, i, epoch)
% returns indeces from vector i which correspond to timestamps in ts within
% epoch start and stop times t_start and t_stop
%
% example: we have indeces i that refer to vector ts of timestamps, give me just the i which
% lie betweeh t_start and t_stop?
%
% also: i_i are indeces in i for which ts is good (used in theta vector

% andrew 18 march 2010

if nargin<3
    error('IsolateEpoch: not enough input arguments');
end

if max(i)>length(ts)
    
    warning('CMBH:error', 'Your indeces do not align with the timestamp vector');
    return
    
end

t = ts(i);

if numel(epoch)>2
    epoch = mat2cell(epoch, ones(1, size(epoch,1)), 2); 

    i_i = cellfun(@(a) find(t>=a(1) & t<=a(2)), epoch, 'Unif', 0);
    
    i_ts = cellfun(@(c) i(c), i_i, 'unif', 0);
    
else
    t_start = epoch(1);
    t_stop = epoch(2);
    
    i_i = find(t>=t_start & t<=t_stop);
    
    i_ts = i(i_i);
    
end
end

function self = SetLFPInds(self)

    if ~iscell(self.p_ind)
        self.p_lfp_ind = cellfun(@(c) nan(size(c)), {self.p_ind}, 'unif', 0); % will eventually be the aligned indices in the LFP
    else
        self.p_lfp_ind = cellfun(@(c) nan(size(c)), self.p_ind, 'unif', 0); % will eventually be the aligned indices in the LFP
    end
                    
end

function self = SetCelInds(self)

    cel = self.cel;
    
    ind = cell(size(self.epoch,1), size(cel,1));
    
    spkind = cell(size(self.epoch,1), size(cel,1));
    
    ind_lfp = cell(size(self.epoch,1), size(cel,1));

    for j = 1:size(cel,1) % now set root.p_cel_ind

        [i, i2] = IsolateEpoch(self.b_ts, self.spike(cel(j,1), cel(j,2)).i, self.epoch);

        if iscell(i)

            if ~isempty(self.active_lfp)
                ind_lfp(:, j) = cellfun(@(c) self.spike(cel(j,1), cel(j,2)).i_lfp{self.active_lfp}(c), i2, 'unif', 0); 
            else
                ind_lfp(:, j) = cell(length(i),1);
            end

            ind(:, j) = i; % video indices
            
            spkind(:,j) = i2;

        else

            if ~isempty(self.active_lfp)
                ind_lfp{1, j} = self.spike(cel(j,1), cel(j,2)).i_lfp{self.active_lfp}(i2); 
            else
                ind_lfp{1, j} = [];
            end

            ind{1, j} = i; % video indices
            
            spkind{1,j} = i2;

        end

    end % now we have cell arrays like {epochs, cells} of indices in video tracking and lfp

    self.p_cel_ind = ind; % cell array like {epochs, cells} of tracking alignment

    self.p_cel_lfp_ind = ind_lfp; % cell array like {epochs, cells} of LFP alignment
            
    self.p_cel_spkind = spkind;
    
end

function l = CheckBaseVarLength(self)

    l(1) = length(self.p_b_x);
    l(2) = length(self.p_b_y);
    l(3) = length(self.p_b_ts);
    l(4) = length(self.p_b_vel);
    l(5) = length(self.p_b_headdir);
    l(6) = length(self.p_b_myvar);

    l = max(l);

end