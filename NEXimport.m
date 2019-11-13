     function varargout = NEXimport(varargin)
            
            % to do: head direction vector?
            
            import CMBHOME.PLX2Mat.* % imports NL import utilities for windows only
            import CMBHOME.Utils.*
            import CMBHOME.Session
            import CMBHOME.Spike

            p = inputParser;

            p.addParamValue('base_path', '', @(c) ischar(c));
            p.addParamValue('fix_pos', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('fix_headdir', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('batch', 0, @(c) numel(c)==1 && (c==1 || c==0));
            p.addParamValue('save_file', 1, @(c) numel(c)==1 && (c==1 || c==0));
           
            p.parse(varargin{:});

            base_path = p.Results.base_path;
            fix_pos = p.Results.fix_pos;
            fix_headdir = p.Results.fix_headdir;
            batch = p.Results.batch;
            save_file = p.Results.save_file;          
            
            [fopen, fsave] = CollectNEXFiles(batch,base_path); % cell array of directories/filenames to be imported
            
            if save_file, ca_roots = cell(size(fopen,1), 1); end
            
            for ii = 1 : size(fopen,1) % make objects, one by one
                                
                % this script reads all the spike timestamp and a/d info from a plx file into matlab
                % variables.

               
                [nexFile] = readNexFile(fopen{ii,:});
                OpenedFileName = fopen{ii,:};
                Version = nexFile.version;
                freq = nexfile.freq;
                comment = nexFile.comment;
              
                
                
                disp(['Opened File Name: ' OpenedFileName]);
                disp(['Version: ' num2str(Version)]);
                disp(['Frequency : ' num2str(Freq)]);
                disp(['Comment : ' Comment]);
             
             
                % some of the information is only filled if the plx file version is >102
                if ( Version > 102 )
                    if ( Trodalness < 2 )
                        disp('Data type : Single Electrode');
                    elseif ( Trodalness == 2 )
                        disp('Data type : Stereotrode');
                    elseif ( Trodalness == 4 )
                        disp('Data type : Tetrode');
                    else
                        disp('Data type : Unknown');
                    end
                end   
                
                % get some counts
             
                evcounts = cellfun(@(a) length(a.timestamps),nexFile.events);
                evcounts = [evcounts;cellfun(@(a) length(a.timestamps),nexFile.markers) ];
                tscounts = zeros(26,24);
                for jj = 1:length(nexFile.neurons)
                    tet = nexFile.neurons{jj}.wireNumber+1;
                    unit = nexFile.neurons{jj}.unitNumber;
                    tscounts(unit,tet) = length(nexFile.neurons{jj}.timestamps);
                    
                end
                    
                    
           

                % tscounts, wfcounts are indexed by (unit+1,channel+1)
                % tscounts(:,ch+1) is the per-unit counts for channel ch
                % sum( tscounts(:,ch+1) ) is the total wfs for channel ch (all units)
                % [nunits, nchannels] = size( tscounts )
                % To get number of nonzero units/channels, use nnz() function

                % get events
                
                % header information, including event strings and
                % corresponding channel numbers
                headerinf = readPLXHeaders(OpenedFileName);                
                
               
                event_channels = [1:length(nexFile.events) 257];
                evnames = cellfun(@(a) a.name,nexFile.events,'uni',0); % channel names
                
                evnames(event_channels>256) = []; % only want those less than 257
                event_channels(event_channels>256) = [];
                
                event_channels(strcmp(evnames, 'Frame Marker')) = []; % get rid of frame markers as well
                evnames (strcmp(evnames, 'Frame Marker')) = []; 
                
                for iev = 1:299
                    if ( evcounts(iev) > 0 )
                        if ( iev == 257 )
                            % treat strobed channel seperately, just to avoid setting up a
                            % cell array to hold strobed values when only one channel will
                            % have them.
                            [nevs{iev}, tsevs{iev}, svStrobed] = plx_event_ts(OpenedFileName, iev); 
                        else
                            [nevs{iev}, tsevs{iev}, svdummy] = plx_event_ts(OpenedFileName, iev);
                        end
                    end                    
                end
                          
                
                
                for iev = 1:length(event_channels)
                    if ( evcounts(iev) > 0 )
                        if ( event_channels(iev) == 257 )
                            % treat strobed channel seperately, just to avoid setting up a
                            % cell array to hold strobed values when only one channel will
                            % have them.
                            [nevs{iev}, tsevs{iev}, svStrobed] = plx_event_ts(OpenedFileName, iev); 
                        else
                            [nevs{iev}, tsevs{iev}, svdummy] = plx_event_ts(OpenedFileName, iev);
                        end
                    end                    
                end
                
                
                x = [];
                y = [];
                t = [];
                headdir = [];
                
                % build event cell array
                event = cell(sum([nevs{event_channels}]), 2);

                if ~isempty(event_channels)
                    tmp_ind = 1;
                    for i_ch = 1:length(event_channels)
                        if nevs{event_channels(i_ch)}>0
                            event(tmp_ind:tmp_ind+nevs{event_channels(i_ch)}-1,1) = evnames(i_ch); % set labels
                            
                            event(tmp_ind:tmp_ind+nevs{event_channels(i_ch)}-1,2) = num2cell(tsevs{event_channels(i_ch)}); % set times
                                                        
                            tmp_ind = tmp_ind+nevs{event_channels(i_ch)}; % update index
                        end
                    end
                    
                    [~, tmpind] = sort([event{:,2}]); % sort by time

                    event = event(tmpind, :);
                    
                else
                    event = cell(1,2);
                end
                
                % get some other info about the spike channels
                [nspk,spk_filters] = plx_chan_filters(OpenedFileName);
                  [nspk,spk_gains] = plx_chan_gains(OpenedFileName);
                [nspk,spk_threshs] = plx_chan_thresholds(OpenedFileName);
                [nspk,spk_names] = plx_chan_names(OpenedFileName);

                spk_names = cellstr(spk_names);
                
                % get the a/d data into a cell array also.
                % This is complicated by channel numbering.
                % The presence/absence of slow analog data can be seen by looking at the
                % evcounts array at indexes 300-363. E.g. the number of samples for
                % analog channel 0 is stored at evcounts(300).
                % Note that analog ch numbering starts at 0, not 1 in the data, but the
                % 'allad' cell array is indexed by ich+1
                numads = 0;
                for ich = 0:63
                    if ( evcounts(300+ich) > 0 )
                        [adfreq, nad, tsad, fnad, allad{ich+1}] = plx_ad(OpenedFileName, ich);
                        numads = numads + 1;
                    end
                end

                if ( numads > 0 )
                    [nad,adfreqs] = plx_adchan_freqs(OpenedFileName);
                    [nad,adgains] = plx_adchan_gains(OpenedFileName);
                    [nad,adnames] = plx_adchan_names(OpenedFileName);

%                     % just for fun, plot the channels with a/d data
% 
%                     [adrows,nActiveADs] = size(allad);
%                     for ich = 1:nActiveADs
%                         if ( size(allad{ich}) > 0 )
%                             subplot(nActiveADs,1,ich); plot(allad{ich});
%                         end
%                     end
                end
                
                if nevs{257}>0 % check that there is position data
                    
                    [n_coords, trash, nVTMode, c] = plx_vt_interpret(tsevs{257}, svStrobed);    % added by andrew to spit out coordinates, etc
                    
                    % check plx_vt_interpret to learn about nVTMode
                    
                    switch nVTMode
                        case {1,2,3,4,5}
                    
                            y = c(:, 3);
                            x = c(:, 2);
                            t = c(:, 1);
                            
                        case {6,7,8} % there are two sets of coordinates, presumably one for each LED
                            % I take the mean of both sets of coordinates,
                            % so long as one of them isnt 0. if it is, then
                            % i just assume the other is right
                            
                            % if they are both zero, it stays zero and then
                            % hopefully the FixPos can fix it        
                                                        
                            c((c(:,2)==0 & c(:,3)==0),2:3) = NaN; % set zeros to NaN so we ignore them in mean
                            c((c(:,4)==0 & c(:,5)==0),4:5) = NaN;
                            
                            x = nanmean([c(:,2), c(:,4)],2);
                            y = nanmean([c(:,3), c(:,5)],2);
                            
                            x(isnan(x)) = 0;
                            y(isnan(y)) = 0;
                            
                            t = c(:,1);
                                                    
                        otherwise % nothing programmed for 3 leds yet
                            disp('Three LED signals detected, only using the first');
                            y = c(:, 3);
                            x = c(:, 2);
                            t = c(:, 1);

                    end
                    
                else
                    disp('Import.PLX: no tracking data found. This import script assumes you ran CinePlex before import.');
                end

                root = Session('name', [fsave{ii,:}], ...
                                'b_ts', t,'b_x',x, 'b_y', y, 'b_headdir',...
                                nan(length(y),1), 'event', event, 'raw_pos', 1, ...
                                'raw_headdir', 1, 'date_created', now, ...
                                'epoch', [t(1) t(end)], 'fs_video', mean(diff(t)).^-1, ...
                                'path_raw_data', fopen{ii}, 'path_lfp', {[fopen{ii,:}]}); % put all behavioral stuff together

                if fix_pos
                    root = root.FixPos();
                end
                
                if fix_headdir
                    root = root.FixDir();
                end                
                                                 
                clear x y t headdir                
                                
                root.spike = Spike(); % initialize spiking data
                
                % gives actual number of units (including unsorted) and actual number of
                % channels plus 1
                [nunits1, nchannels1] = size( tscounts );

                root.spike(1:(size(tscounts,2)-1)/Trodalness, 1:nunits1-1) = Spike(); % initialize the spike array

                % we will read in the timestamps of all units,channels into a two-dim cell
                % array named allts, with each cell containing the timestamps for a unit,channel.
                % Note that allts second dim is indexed by the 1-based channel number.
                for iunit = 1:nunits1-1   % starting with unit 1 (sorted). 0 is unsorted
                    ich = 1;
                    while ich<=nchannels1-1
                        if ( tscounts( iunit+1 , ich+1 ) > 0 )
                            % get the timestamps for this channel and unit
                            %[nts, allts{iunit+1,ich}] = plx_ts(OpenedFileName, ich , iunit );
                            
                            if lower(spk_names{ich}(1))=='t' % this is a tetrode channel
                                tet_ind = sscanf(spk_names{ich}, '%1s%f%1s%f');
                                tet_ind = tet_ind(2);
                                
                                [trash, spk_ts] = plx_ts(OpenedFileName, ich , iunit );
                                
                                spk_i = shiftdim(SpkInds(root.b_ts, spk_ts));
                                                               
                                if isempty(root.spike(tet_ind, iunit).ts) && ~isempty(spk_i) % if this cell hasnt been loaded (because this cells spikes are repeated four times over in tetrode event land)

                                    disp([spk_names{ich} ' is tetrode ' int2str(tet_ind) ' and cell ' int2str(iunit)]);
                                    root.spike(tet_ind, iunit) = Spike('i', spk_i, 'ts', spk_ts, 'tet_name', spk_names{ich});
                                    
                                end
                                
                                if isempty(spk_i), disp(['Skipping tetrode ' int2str(tet_ind) ' and cell ' int2str(iunit) ' because there were no spikes.']); end
                                
                                ich = ich+1; % skip next three channels, because they are just the identical spikes from other electrodes in the tetrode
                            else
                                ich = ich+1;
                            end
                        else
                            ich = ich+1;
                        end
                    end
                end
                
                if save_file
                    save([fsave{ii,:}], 'root'); % save
                end
                
                if nargout==1
                    ca_roots{ii} = root;
                end
                
            end           
        
            if nargout==1
                varargout{1} = ca_roots;
            end
            
            function [fopen, fsave] = CollectNEXFiles(batch, base_path)

                fopen = cell(1,2); % reasonably large cell array of str fnames
                fsave = cell(1,2);

                if batch        % if batch is selected, process each folder in base_path
                                % if directory cluster_files exists in each folder,
                                % no user input will be necessary
                    if isempty(base_path), error('For batch mode you must pass the base_path argument'); end

                    folders = dir(base_path);
                    
                    folders(end+1).name = base_path; % check base_path, too

                    f_ind = 1;

                    for folderi = 1:length(folders)
                        if folders(folderi).isdir & ~strcmp(folders(folderi).name,'.') & ~strcmp(folders(folderi).name,'..') % for every folder, call Import.NL

                            tmp_base_path = fullfile(base_path, folders(folderi).name);

                            tmp_file = dir([tmp_base_path '*.nex']);
                            
                            tmp_file = {tmp_file(:).name};
                            
                            for i=1:length(tmp_file)

                                fopen(i, :) = {tmp_base_path, tmp_file{i}};
                                fsave(i,:) = {tmp_base_path, ['CMBH_', strrep(tmp_file{i}, '.nex', '.mat')]};
                                
                            end
                       
                        end
                    end
                                        
                    return;
                end % end batch 

                % otherwise, open base_path (pwd if empty) and keep
                % returning files until user says done

                if isempty(base_path), base_path = pwd; end

                f_ind = 1;

                while 1 % loop until user cancels, and we get no files

                    [load_files, base_path] = uigetfile('*.plx','Please select PLX files to load. Exit to finish.',base_path, 'MultiSelect', 'on');

                    if iscell(load_files)
                            for i = 1:length(load_files)
                                fopen(f_ind,:) = {base_path, load_files{i}};
                                fsave(f_ind,:) = {base_path, ['CMBObj_' strrep(load_files{i}, '.plx', '.mat')]};
                                f_ind = f_ind+1;
                            end
                    elseif ischar(load_files)
                            fopen(f_ind,:) = {base_path, load_files};
                            fsave(f_ind,:) = {base_path, ['CMBObj_' strrep(load_files, '.plx', '.mat')]};
                            f_ind = f_ind+1;
                    else
                        break; % exit loop collecting files
                    end

                end

                if isempty(fopen)
                    disp('No .plx files found/selected. No CMBObjects will be saved.');
                end

            end   
        end