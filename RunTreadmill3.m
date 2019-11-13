function RunTreadmill3(settings)
% $Id: RunTreadmill.m 4831 2013-05-13 19:56:09Z nrobinson $

rstream = RandStream('mt19937ar','Seed', mod(prod(clock()),2^32));
RandStream.setGlobalStream(rstream);

tr= Treadmill('COM4');
laserserial = serial('COM3','BaudRate',9600);
fopen(laserserial);
autocloselaserserial = onCleanup(@() fclose(laserserial));

defsettings.mapserver = false;
defsettings.joystick = true;
defsettings.tracking = false;
defsettings.maze = 'Treadmill';
defsettings.training= false;
defsettings.deltaspeed = 5;
defsettings.laserdelay = 2; % seconds till laser on
defsettings.laserduration = 2; % seconds
defsettings.lasermode = 'nolaser'; % 'nolaser' = no laser
                                   % '5050'    = 50% no laser, 50% laser
                                   % 'arm'     = 20% arm, 40% laser, 40% no laser

defsettings.treadmill.on = true;
defsettings.treadmill.speeds = 30;
defsettings.treadmill.times = 10 ;
defsettings.treadmill.distances = 1500;
defsettings.treadmill.fixed = 'time'; % Fixed parameter
defsettings.treadmill.vary = 'speed'; % Independant parameter

defsettings.treadmill.pausenotstop = true; % Remember settings when manually stopping treadmill

tim = timer('StartFcn',@startfcn,'TimerFcn',@timerfcn,...
    'StopFcn',@stopfcn,'ErrorFcn',@errorfcn,'Period',0.01,...
    'ExecutionMode','fixedDelay','Name','Main RunTreadmill Timer');

s = 0;
deviceNum = 0;
needtosave = false;

EVENTCODES = maze_events();       % Load the event codes.

history = [];

timeref = tic;

lap=1;
sessionnum = 0;
session.comment = '';
session.date = now();
session.laps = [];
session.correct = [];
session.dig = [];
session.objects = [];
session.laser = [];
session.samples = [];

trialblocksize = 16;
trialtypes = [];
objectlist = [];
laserlist = [];

treadmillstart = 0;
treadmillcount = 0;
treadmillind = 0;
treadmillcurspeed = 0;
treadmillcurtime = 0;
treadmillpausestate = 0;

trackstat.ready = false;
joystatus = 0;

if(nargin == 0)
    settings = checksettings();
else
    settings = checksettings(settings);
end
savetrackstat = trackstat; % DEBUGGING

h = waterportsGUI();
set(h.controls,'CloseRequestFcn',@closeButton)
set(h.timer,'CloseRequestFcn',@toggleTimer);
% set(h.c,'Callback',{@(x,y) process_event(EVENTCODES.,'computer')});
% set(h.i,'Callback',{@(x,y) process_event(EVENTCODES.,'computer')});
set(h.tmillstart,'Callback',{@(x,y) process_event(EVENTCODES.TM,'computer')});
set(h.tmillstop,'Callback',{@(x,y) process_event(EVENTCODES.TMS,'computer')});
set(h.open,'Callback',@openHist);
set(h.save,'Callback',@saveHist);
set(h.pause,'Callback',@playpause);
set(h.reset,'Callback',@resetchecker);
set(h.settings,'Callback',@settingscallback);
set(h.showtimer,'Callback',@toggleTimer);
applysettings;

    function startfcn(varargin)
        output(h.maindsp, 'Disabling Screen Saver...');
        [status, result] = dos('FlipSS /off');
        if(status == 0); output(h.maindsp, 'Screen Saver Disabled.');
        else output(h.maindsp, regexprep(result,'\n',' '));
        end
        
        % Open connection with the Plexon Server
        if(settings.mapserver)
            output(h.maindsp,'Connecting to Plexon MAP Server...');
            try s = PL_InitClient(0);
            catch %#ok<CTCH>
                s = 0;
            end
            
            if s == 0
                output(h.maindsp,'Failed to connect to Plexon MAP Server');
            end
        end
        
        if(s ~= 0)
            output(h.maindsp, 'Connected to Plexon MAP Server.');
            
            pars = PL_GetPars(s);
            if(~pars(11))
                output(h.maindsp, 'Sort Client is not running.');
                warndlg('Sort Client must be running to use the button box or video tracking. Disconnecting from MAP Server.','Sort Client not running.');
                output(h.maindsp, 'Closing connection to Plexon MAP Server...');
                PL_Close(s);
                s = 0;
            end
        end
        
        set(h.pause,'String','Pause');
        set(h.pause,'Enable','on');
    end

    function errorfcn(varargin)
        output(h.maindsp, 'Error encountered.');
    end

    function stopfcn(varargin)
        output(h.maindsp, sprintf('Pausing...'));
        
        if (s ~= 0);
            output(h.maindsp, 'Closing connection to Plexon MAP Server...');
            PL_Close(s);
            s = 0;
        end
        
        output(h.maindsp, 'Enabling Screen Saver...');
        [status, result] = dos('FlipSS /on');
        if(status == 0); output(h.maindsp, 'Screen Saver Enabled.');
        else output(h.maindsp, regexprep(result,'\n',' '));
        end
        
        set(h.pause,'String','Start');
        set(h.pause,'Enable','on');
    end

    function timerfcn(varargin)
        if(settings.joystick && isstruct(joystatus) )
            joyevents = gamepad('poll');
            for k = 1:size(joyevents,1)
                joystick_event(joyevents(k,1),joyevents(k,2));
            end
        end
        
        if(s ~= 0)
            res = PL_WaitForServer(s, 10);
            if(res ~= 0)
                [~, t] = PL_GetTS(s);
                
                if(settings.tracking)
                    coords = PL_GetCoords(t);
                    
                    savetrackstat = trackstat; % DEBUGGING
                    
                    [trackstat,trackevents] = tracking('track',trackstat,coords);
                    
                    % DEBUGGING
                    if(size(trackevents,1)> 0); disp(trackevents); end
                    if(isfield(trackstat,'history') && isfield(savetrackstat,'history'))
                        if(any(trackstat.history ~= savetrackstat.history))
                            if(length(trackstat.history) >= 10)
                                disp(trackstat.history(end:-1:end-9));
                            else
                                disp(trackstat.history(end:-1:1));
                            end
                        end
                    end
                    % END DEBUGGING
                    
                    for j = 1:size(trackevents,1);
                        process_event(trackevents(j,2),'tracking',now());
                    end
                end
            end
        end
        updateTimers();
    end

    function joystick_event(event_num,event_time)
        jsmap = joystickmap;
        ind = find(event_num == jsmap(:,1));
        if(~isempty(ind))
            event_num = jsmap(ind,2);
            process_event(event_num,'button',event_time);
        end
    end

    function process_event(event_num,~,event_time)
        if(nargin~=3); event_time = now(); end
        sessionts = (event_time - session.date)*24*60*60;
        
        switch event_num
            case EVENTCODES.TM
                output(h.maindsp, sprintf('Starting treadmill (%.2f)', sessionts));
                initTreadmill(event_num,event_time);
                
                if(s ~= 0); PL_SendUserEvent(s, event_num); end
            case EVENTCODES.TMS
                output(h.maindsp, sprintf('Stopping treadmill (%.2f)', sessionts));
                stopTreadmill();
                if(s ~= 0); PL_SendUserEvent(s, event_num); end
            case EVENTCODES.DECSPEED
                if(settings.training~=0)
                    settings.treadmill.speeds = settings.treadmill.speeds -settings.deltaspeed;
                    if(settings.treadmill.speeds < 0); settings.treadmill.speeds = 0; end
                    output(h.maindsp, sprintf('Speed= %d', settings.treadmill.speeds));
                    tr.SetSpeed=settings.treadmill.speeds;
                end
            case EVENTCODES.INCSPEED
                if(settings.training~=0)
                    settings.treadmill.speeds = settings.treadmill.speeds +settings.deltaspeed;
                    if(settings.treadmill.speeds > tr.MaxSpeed); settings.treadmill.speeds = tr.MaxSpeed; end
                    output(h.maindsp, sprintf('Speed= %d', settings.treadmill.speeds));
                    tr.SetSpeed=settings.treadmill.speeds;
                end
            case EVENTCODES.SAMPLE
                if(numel(session.samples)==0);
                    needtosave = true;
                    timeref = tic;
                    session.date = event_time;
                    sessionts = (event_time - session.date)*24*60*60;
                end
                session.samples(end+1) = event_time;
                output(h.maindsp, sprintf('Sample (%.2f)', sessionts));
                if(s ~= 0); PL_SendUserEvent(s, event_num); end
            case {EVENTCODES.DIG, EVENTCODES.NODIG}
                ratresponse(event_num, event_time)
                if(s ~= 0); PL_SendUserEvent(s, event_num); end
            case EVENTCODES.LASERZONE
                armlaser(event_num, event_time)
                if(s ~= 0); PL_SendUserEvent(s, event_num); end
        end
        updateCounter();
        updateDisplay();
    end

    function createtrialtypes()
        if(settings.laserduration <= 0 && strcmpi(settings.lasermode,'5050'))
            settings.lasermode = 'nolaser';
        end

        switch settings.lasermode
            case 'arm'
                trialtypes = mod(randperm(trialblocksize),20);
                objectlist = (mod(trialtypes,2)==0);
                
                laserlist = nan(size(trialtypes));
                laserlist(trialtypes<8) = 0;                  % no laser
                laserlist(trialtypes>=8 & trialtypes<16) = 1; % treadmill laser
                laserlist(trialtypes>=16) = 2;                % arm laser
            case '5050'
                trialtypes = mod(randperm(trialblocksize),8);
                objectlist = (mod(trialtypes,2)==0);
                laserlist = (trialtypes>3);
            otherwise
                trialtypes = mod(randperm(trialblocksize),4);
                objectlist = (mod(trialtypes,2)==0);
                laserlist = false(size(trialtypes));
        end
        
        [objectlist]=trialGenCM(objectlist');
    end
        
    function ratresponse(event_num, event_time)
        needtosave = true;
        lap_ = mod(lap-1,numel(objectlist))+1;
        object = objectlist(lap_);
        laser = laserlist(lap_);
        
        session.laps(lap) = event_time;
        session.objects(lap) = object;
        session.laser(lap) = laser;
        
        sessionts = (event_time - session.date)*24*60*60;
        
        % Determine whether he should have dug or not:
        if(object == 0); shoulddig = true;
           
        else
            
            if(object == 1); shoulddig = false;
            end  
        end
            
        switch event_num
            case EVENTCODES.DIG
                session.dig(lap) = true;
                if(shoulddig)
                    session.correct(lap) = true;
                    output(h.maindsp, sprintf('Correct Dig (%.2f)',sessionts));
                else
                    session.correct(lap) = false;
                    output(h.maindsp, sprintf('Incorrect Dig (%.2f)',sessionts));
                end
            case EVENTCODES.NODIG
                session.dig(lap) = false;
                if(shoulddig)
                    session.correct(lap) = false;
                    output(h.maindsp, sprintf('Incorrect No Dig (%.2f)',sessionts));
                else
                    session.correct(lap) = true;
                    output(h.maindsp, sprintf('Correct No Dig (%.2f)',sessionts));
                end
        end
                
        lap = lap+1;
        
        if(mod(lap-1,numel(objectlist)) == 0)
            createtrialtypes;
        end
    end
        
    function armlaser(event_num, event_time)
        lap_ = mod(lap-1,numel(objectlist))+1;
        if(~isempty(laserlist) && lap > 0 && laserlist(lap_)==2 &&...
                settings.laserduration > 0)
            pulselaser(laserserial,0,settings.laserduration);
            output(h.maindsp, 'Pulsing Laser (arm)');
        end
    end

        % initTreadmill starts the treadmill sequence, which could include a
        % waterport opening, followed by a pause, followed by the actual start
        % of the treadmill.
        function initTreadmill(event_num, event_time)
            if(strcmp(tr.Running,'on'))
                output(h.maindsp, 'Treadmill is already running.');
            else
                needtosave = true;
                tr.SetSpeed=settings.treadmill.speeds;
                tr.start(settings.treadmill.times)
                lap_ = mod(lap-1,numel(objectlist))+1;
                if(~isempty(laserlist) && lap > 0 && laserlist(lap_)==1 &&...
                        settings.laserduration > 0)
                    pulselaser(laserserial,settings.laserdelay,settings.laserduration);
                    output(h.maindsp, 'Pulsing Laser (treadmill)');
                end
            end
        end
        
        % stopTreadmill is called to force the treadmill to stop immediately.
        function stopTreadmill()
            if(strcmp(tr.Running,'on'))
                if(settings.treadmill.pausenotstop)
                    treadmillpausestate = 1;
                    output(h.maindsp, 'Treadmill sequence forced to pause');
                else
                    output(h.maindsp, 'Treadmill sequence forced to stop');
                end
                if(settings.treadmill.pausenotstop)
                    treadmillpausestate = 2;
                    output(h.maindsp, 'Treadmill paused');
                end
            end
            tr.stop();
        end
        
        function output(display, str)
            output_list = get(display,'String');
            output_list = output_list(end:-1:1);
            num = length(output_list)+1;
            str = strtrim(str);
            output_list{end+1} = sprintf('%3d  %s',num,str);
            set(display,'String',output_list(end:-1:1));
            drawnow
        end
        
        function settingscallback(varargin) %#ok<VANUS>
            settings = changesettings(settings, h);
            settings = checksettings(settings);
            applysettings;
        end
        
        function applysettings()
            if(strcmp('treadmill',regexprep(lower(settings.maze),'[^a-z]','')))
                set(h.tmillstart,'Enable','on');
                set(h.tmillstop,'Enable','on');
            else
                set(h.tmillstart,'Enable','off');
                set(h.tmillstop,'Enable','off');
            end
            
            if(isempty(laserlist) || ...
                    (any(laserlist==2) && ~strcmpi(settings.lasermode,'arm')) || ...
                    (all(laserlist==0) && ~strcmpi(settings.lasermode,'nolaser')) || ...
                    (any(laserlist==1)&& all(laserlist<2) && ~strcmpi(settings.lasermode,'5050')) || ...
                    (all(~laserlist) && settings.laserduration>0) || ...
                    (any(laserlist) && settings.laserduration==0))
                createtrialtypes;
            end
            
            updateDisplay();
            updateCounter();
        end
        
        function settings = checksettings(settings)
            if(nargin==0); settings = defsettings;
            else
                % Now check all the rest of the settings, making sure they are
                % at least present and the correct variable type.
                settings = checkfieldtypes(settings, defsettings);
            end
            
            % If Joystick control is activated, then initialize the joystick
            % listener. Otherwise stop it and delete it.
            if(settings.joystick)
                if(~isstruct(joystatus));
                    try
                        list = gamepad('list');     % List IDs of plugged in gamepads
                       gamepad('start',list);        % Start polling gamepad ID
                        joystatus = gamepad('status'); % Output the gamdpad status structure
                        gamepad('stop');            % Stop polling gamepad
                        events = gamepad('poll');   % Get list of events in the queue
                        
                    catch %#ok<CTCH>
                        status = 0;
                        warndlg('Joystick failed to initialize.','Joystick Disabled');
                        settings.joystick = false;
                        
                    end
                end
         
            end
            
            % If MAP Server access is disabled, then disable tracking.
            if(~settings.mapserver && settings.tracking)
                warndlg('Tracking will not work with MAP Server access disabled.','Tracking Disabled');
                settings.tracking = false;
            end
            
            % If tracking is currently off, then reset the tracking
            % information, leave trackstat.ready = false.
            if(~settings.tracking)
                tmptrackstat.ready = false;
                trackstat = tmptrackstat;
                clear tmptrackstat

            % If tracking is currently on, but trackstat.ready isn't a field,
            % then re-initalize the tracking status.
            elseif(~isfield(trackstat,'ready') || ~islogical(trackstat.ready))
                trackstat = tracking('initialize',settings.maze);

            % If tracking is currently on, and the trackstat variable has valid
            % data, then keep that data, but validate it.
            else
                trackstat = tracking('validate',trackstat,settings.maze);
            end

            % If trackstat.ready is false, then disable tracking.
            if(settings.tracking && ~trackstat.ready);
                settings.tracking = false;
                warndlg('Tracking could not initialize.','Tracking Disabled');
            end

            % Check the treadmill settings
            switch settings.treadmill.fixed
                case 'time'
                    settings.treadmill.times = settings.treadmill.times(1);
                    % Time Fixed, Speed Varying, Computer Time
                    settings.treadmill.vary = 'speed';
                    settings.treadmill.speeds = sort(settings.treadmill.speeds);
                    settings.treadmill.distances = settings.treadmill.speeds.*settings.treadmill.times;
                    
                    
                    settings.treadmill.fixed = defsettings.treadmill.fixed;
                    
                    settings.treadmill.times = settings.treadmill.distances./settings.treadmill.speeds;
            end
            
            if(max([length(settings.treadmill.times) length(settings.treadmill.speeds)])<treadmillind)
                treadmillcount = 0;
                treadmillind = 0;
            end
        end
        
        function newstruct = checkfieldtypes(newstruct, refstruct)
            fields = fieldnames(refstruct);
            
            for j = 1:size(fields,1)
                if(~isfield(newstruct,fields{j}) || ...
                        ~isa(newstruct.(fields{j}),class(refstruct.(fields{j}))))
                    newstruct.(fields{j}) = refstruct.(fields{j});
                elseif(isnumeric(newstruct.(fields{j})) && ...
                        (isempty(newstruct.(fields{j})) || ...
                        any(isnan(newstruct.(fields{j})))))
                    newstruct.(fields{j}) = refstruct.(fields{j});
                    warning(['Invalid number entered for ' fields{j} ', resetting to default value.'],'Invalid Entry');
                elseif(isstruct(newstruct.(fields{j})))
                    newstruct.(fields{j}) = checkfieldtypes(newstruct.(fields{j}), refstruct.(fields{j}));
                end
            end
        end
        
        function openHist(varargin) %#ok<VANUS>
            if(~needtosave || (needtosave && resetchecker(true)))
                historytmp = loadHistory();
                if(~isempty(historytmp))
                    history = historytmp;
                    sessionnum = 0;
                end
                
                if(~isempty(history))
                    if(isfield(history,'sessions') && ...
                            isfield(history.sessions,'settings') && ...
                            ~isempty(history.sessions(end).settings))
                        settings = checksettings(history.sessions(end).settings);
                    elseif(isfield(history,'settings') && ...
                            ~isempty(history.settings))
                        settings = checksettings(history.settings);
                    end
                    output(h.maindsp, sprintf('History for %s loaded.',history.ratid));
                    applysettings;
                end
            end
        end
        
        function success = saveHist(varargin) %#ok<VANUS>
            session.TreadmillLog = tr.Log;
            
            if(needtosave)
                if(~isfield(session,'comment'))
                    session.comment = char(inputdlg({'Session Comments'},...
                        'Session Comments',10));
                else
                    session.comment = char(inputdlg({'Session Comments'},...
                        'Session Comments',10,{session.comment}));
                end
                
                session.maze = settings.maze;
                session.settings = settings;
                
                [history, success, sessionnum] = addSession(session, history, sessionnum);
                
                if(success)
                    [success, history] = saveHistory(history);
                    if(success)
                        output(h.maindsp, sprintf('History for %s saved.',history.ratid));
                        needtosave = false;
                    end
                else
                    sessionnum = 0;
                end
                
            else msgbox('Nothing to save.');
            end
            updateDisplay();
        end
        
        function updateDisplay()
            settingsstr = sprintf('');
            
            if(settings.treadmill.on); settingsstr = [settingsstr ', Treadmill: on'];
            else settingsstr = [settingsstr ', Treadmill: off'];
            end
            
            set(h.settingsdsp,'String',settingsstr);
            
            if(lap < 0)
                sessionstr = sprintf('Not Started');
            elseif(lap == 0)
               % sessionstr = sprintf('Lap: %d, T1: %g, T2: %g',lap,tim.AveragePeriod,joy.AveragePeriod);
            else
                sessionstr = 'Lap: %d, Correct: %d (%.2f%%), T1: %g, T2: %g';
                %sessionstr = sprintf(sessionstr,lap,sum(session.correct),...
                 %   sum(session.correct)/lap*100,tim.AveragePeriod,joy.AveragePeriod);
            end
            
            if(~isempty(history))
                sessionstr = ['Rat ID: ' history.ratid ', ' sessionstr];
            end
            
            set(h.sessiondsp,'String',sessionstr);
            
            if(needtosave)
                set(h.save,'Enable','on');
            else
                set(h.save,'Enable','off');
            end
        end
        
        function updateTimers()
            ts = toc(timeref);
            if(ishandle(h.timerstr))
                timestr = sprintf('%2.0f:%02.0f',floor(ts/60),mod(floor(ts),60));
                set(h.timerstr,'String',timestr);
            end
        end
        
        function updateCounter()
            if(ishandle(h.counter))
                % Automatically loop back to the beginning of the list once
                % you reach the end.
                lap_ = mod(lap-1,numel(objectlist))+1;
                
                if(laserlist(lap_)==1); lstr = 'L';
                elseif(laserlist(lap_)==2); lstr = 'LA';    
                else lstr = ' ';
                end
                
                % Display lap number, object
                set(h.counter,'String',sprintf('%3d: %1d %1s %2s',lap,objectlist(lap_),lstr));
            end
        end
        
        function success = resetchecker(varargin)
            if(islogical(varargin{1}) && varargin{1})
                button = 'Yes';
            else
                button = questdlg('Are you sure you want to reset?','Reset Confirmation...','Yes','No','Yes');
            end
            
            if(strcmp(button,'Yes') && promptSaveSession)
                success = reset();
            else
                output(h.maindsp, 'Save failed or reset cancelled');
                success = false;
            end
        end
        
        function success = reset()
            needtosave = false;
            
            try
                tr.stop();
            catch ME
                fprintf(2,'Error stopping treadmill during reset.\n');
                fwrite(1,ME.getReport());
            end
            tr.clearLog();
            
            timeref = tic;
            
            lap = 1;
            sessionnum = 0;
            session.comment = '';
            session.date = now();
            session.laps = [];
            session.correct = [];
            session.dig = [];
            session.objects = [];
            session.laser = [];
            session.samples = [];
            
            trialtypes = [];
            objectlist = [];
            laserlist = [];
            createtrialtypes;
            
            treadmillstart = 0;
            treadmillcount = 0;
            treadmillind = 0;
            treadmillcurspeed = 0;
            treadmillcurtime = 0;
            treadmillpausestate = 0;
            
            % If tracking is on, we want to clear the tracking information but
            % leave trackstat ready for tracking.
            if(settings.tracking)
                trackstat = tracking('initialize',settings.maze);
                
                % If tracking is off, we want to clear the tracking information,
                % but we don't need to bother trying to reinitialize trackstat.
            else
                tmptrackstat.ready = false;
                trackstat = tmptrackstat;
                clear tmptrackstat
            end
            
            % If trackstat.ready is false, then disable tracking.
            if(settings.tracking && ~trackstat.ready);
                settings.tracking = false;
                warndlg('Tracking could not reinitialize.','Tracking Disabled');
            end
            
            set(h.maindsp,'String',{'  1  Session Reset.'});
            
            if(ishandle(h.counter))
                set(h.timerstr,'String','');
                set(h.counter,'String','');
                set(h.timerstr,'ForegroundColor','white');
                set(h.counter,'ForegroundColor','white');
            end
            
            updateDisplay();
            updateCounter();
            success = true;
        end
        
        function success = promptSaveSession()
            if(needtosave);
                button = questdlg('Save unsaved session?','Saving Session...');
                switch button
                    case {'Yes'}
                        success = saveHist();
                    case {'No'}
                        button2 = questdlg('Are you sure? You will lose your unsaved session.','Clear Session Confirmation...','Yes','No','No');
                        switch button2
                            case {'Yes'}
                                success = true;
                            otherwise
                                success = false;
                        end
                    otherwise
                        success = false;
                end
            else
                success = true;
            end
        end
        
        function playpause(varargin) %#ok<VANUS>
            % If currently running, then stop.
            if(strcmp(tim.Running,'on'))
                set(h.pause,'Enable','off');
                stop(tim);
                % If stopped, then start running.
            else
                set(h.pause,'Enable','off');
                start(tim);
            end
        end
        
        function closeButton(source, varargin) %#ok<VANUS>
            if(strcmp(tim.Running,'off') && promptSaveSession)
              
                delete(source);
                if(ishandle(h.timer))
                    delete(h.timer);
                end
                if(isvalid(tim))
                    delete(tim);
                end
                
                try
                    tr.delete;
                catch ME
                    fprintf(2,'Manually closing serial port.\n');
                    fwrite(1,ME.getReport());
                    Treadmill.deleteSerial();
                end
            else warndlg('You must click ''Pause'' before closing this window.','Cannot Close');
            end
        end
        
        function toggleTimer(varargin) %#ok<VANUS>
            switch get(h.timer,'Visible');
                case 'on'
                    set(h.timer,'Visible','off');
                    set(h.showtimer,'String','Show Timer');
                case 'off'
                    set(h.timer,'Visible','on');
                    set(h.showtimer,'String','Hide Timer');
                    maximize(h.timer);
            end
        end
        
        function h = waterportsGUI()
            screen = get(0,'ScreenSize');
            winsize = [450 600];
            
            h.controls = figure('Name','Treadmill Task',...
                'NumberTitle','off','MenuBar','none',...
                'Position', [10 (screen(4)-winsize(2)-30) winsize]);
            
            h.buttons = uipanel('Parent',h.controls,...
                'Units','normalized','Position',[0.01 0.82 0.35 0.17]);
            
            
%             h.c = uicontrol(h.buttons,'Style','pushbutton',...
%                 'Units','normalized','Position',[0.02 0.35 0.47 0.30],...
%                 'String','Correct');
%             h.i = uicontrol(h.buttons,'Style','pushbutton',...
%                 'Units','normalized','Position',[0.51 0.35 0.47 0.30],...
%                 'String','Incorrect');
            
            h.open = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.37 0.94 0.20 0.05],...
                'String','Open');
            
            h.save = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.37 0.88 0.20 0.05],...
                'String','Save','Enable','off');
            
            h.tmillstart = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.37 0.82 0.20 0.05],...
                'String','Start Treadmill','Enable','off');
            
            h.pause = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.58 0.94 0.20 0.05],...
                'String','Start');
            
            h.reset = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.58 0.88 0.20 0.05],...
                'String','Reset');
            
            h.tmillstop = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.58 0.82 0.20 0.05],...
                'String','Stop Treadmill','Enable','off');
            
            h.settings = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.79 0.94 0.20 0.05],...
                'String','Settings');
            
            h.showtimer = uicontrol(h.controls,'Style','pushbutton',...
                'Units','normalized','Position',[0.79 0.88 0.20 0.05],...
                'String','Show Timer');
            
            h.settingsdsp = uicontrol(h.controls,'Style','text',...
                'Units','normalized','Position',[0.01 0.79 0.98 0.02]);
            
            h.sessiondsp = uicontrol(h.controls,'Style','text',...
                'Units','normalized','Position',[0.01 0.76 0.98 0.02]);
            
            h.maindsp = uicontrol(h.controls,'Style','edit',...
                'Units','normalized','Position',[0.01 0.01 0.98 0.74],...
                'BackgroundColor','white','HorizontalAlignment','left',...
                'Max',2,'Enable','inactive');
            
            h.timer = figure('Name','Lap Counter',...
                'NumberTitle','off','MenuBar','none','Visible','off',...
                'Position', [5 39 screen(3)-11 screen(4)-63]);
            
            h.counter = uicontrol(h.timer,'Style','text',...
                'Units','normalized','Position',[0 0.0 1 .5],...
                'BackgroundColor','black','ForegroundColor','white',...
                'FontUnits','normalized','FontSize',.65);
            
            h.timerstr = uicontrol(h.timer,'Style','text',...
                'Units','normalized','Position',[0 0.5 1 .5],...
                'BackgroundColor','black','ForegroundColor','white',...
                'FontUnits','normalized','FontSize',.95);
        end
        
        function settings = changesettings(settings, h)
            LASERMODES = {'nolaser','5050','arm'};
            
            currentmode = 1;
            for j = 1:length(LASERMODES)
                if(strcmpi(settings.lasermode,LASERMODES{j})); currentmode = j; break; end
            end
            
            parentwin = get(h.controls,'Position');
            
            setwin = dialog('Name','Settings',...
                'NumberTitle','off','MenuBar','none');
            
            labels.laserdelay = uicontrol(setwin,'Style','text',...
                'String','Laser Delay (s):');
            
            controls.laserdelay = uicontrol(setwin,'Style','edit',...
                'BackgroundColor','white','String',num2str(settings.laserdelay));
            
            labels.laserduration = uicontrol(setwin,'Style','text',...
                'String','Laser Duration (s):');
            
            controls.laserduration = uicontrol(setwin,'Style','edit',...
                'BackgroundColor','white','String',num2str(settings.laserdelay));
            
            labels.lasermode = uicontrol(setwin,'Style','text',...
                'String','Laser Mode:');

            controls.lasermode = uicontrol(setwin,'Style','popupmenu',...
                'String',LASERMODES,'Value',currentmode);

            labels.mapserver = uicontrol(setwin,'Style','text',...
                'String','MAP Server:');
            
            controls.mapserver = uicontrol(setwin,'Style','checkbox','Min',0,'Max',1,...
                'Value',settings.mapserver);
            
            labels.joystick = uicontrol(setwin,'Style','text',...
                'String','Joystick:');
            
            controls.joystick = uicontrol(setwin,'Style','checkbox','Min',0','Max',1,...
                'Value',settings.joystick);
            
            labels.tracking = uicontrol(setwin,'Style','text',...
                'String','Video Tracking:');

            controls.tracking = uicontrol(setwin,'Style','checkbox','Min',0,'Max',1,...
                'Value',settings.tracking);

            labels.treadmillon = uicontrol(setwin,'Style','text',...
                'String','Treadmill:');
            
            controls.treadmillon = uicontrol(setwin,'Style','checkbox','Min',0,'Max',1,...
                'Value',settings.treadmill.on);
            
            labels.tmillspeeds = uicontrol(setwin,'Style','text',...
                'String','Speeds (cm/s):');
            
            controls.tmillspeeds = uicontrol(setwin,'Style','edit',...
                'BackgroundColor','white','String',num2str(settings.treadmill.speeds));
            
            labels.tmilltimes = uicontrol(setwin,'Style','text',...
                'String','Times (s):');
            
            controls.tmilltimes = uicontrol(setwin,'Style','edit',...
                'BackgroundColor','white','String',num2str(settings.treadmill.times));
            
            labels.pausenotstop = uicontrol(setwin,'Style','text',...
                'String','Pause (not Stop):');
            
            controls.pausenotstop = uicontrol(setwin,'Style','checkbox','Min',0,'Max',1,...
                'Value',settings.treadmill.pausenotstop);
            
            labels.training = uicontrol(setwin,'Style','text',...
                'String','Training Mode');
            
            controls.training = uicontrol(setwin,'Style','checkbox','Min',0,'Max',1,...
                'Value',settings.training);
            
            labels.deltaspeed = uicontrol(setwin,'Style','text',...
                'String','Change Speed amount');
            
            controls.deltaspeed = uicontrol(setwin,'Style','edit',...
                'BackgroundColor','white','String',num2str(settings.deltaspeed));
            
            settingnames = fieldnames(controls);
            
            numrows = length(settingnames);
            
            winsize = [200 30+25*numrows];
            set(setwin,'Position',[parentwin(1)+50,...
                parentwin(2)+parentwin(4)-winsize(2)-150,winsize]);
            
            for j = 1:numrows
                ypos = winsize(2)-j*25;
                set(labels.(settingnames{j}),'Position',[05 ypos 90 16]);
                set(controls.(settingnames{j}),'Position',[100 ypos 95 20]);
            end
            
            uicontrol(setwin,'Style','pushbutton','String','OK',...
                'Position',[110 05 30 20],'Callback',{@(x,y) uiresume(setwin)});
            uicontrol(setwin,'Style','pushbutton','String','Cancel',...
                'Position',[145 05 50 20],'Callback',{@(x,y) delete(setwin)});
            
            uiwait(setwin);
            
            if(ishandle(setwin))
                settings.laserdelay = str2num(get(controls.laserdelay,'String')); %#ok<ST2NM>
                settings.laserduration = str2num(get(controls.laserduration,'String')); %#ok<ST2NM>
                settings.lasermode = LASERMODES{get(controls.lasermode,'Value')};
                settings.mapserver = logical(get(controls.mapserver,'Value'));
                settings.joystick = logical(get(controls.joystick,'Value'));
                settings.tracking = logical(get(controls.tracking,'Value'));
                settings.treadmill.on = logical(get(controls.treadmillon,'Value'));
                settings.treadmill.speeds = str2num(get(controls.tmillspeeds,'String')); %#ok<ST2NM>
                settings.treadmill.times = str2num(get(controls.tmilltimes,'String')); %#ok<ST2NM>
                settings.treadmill.pausenotstop = logical(get(controls.pausenotstop,'Value'));
                settings.training = logical(get(controls.training, 'Value'));
                settings.deltaspeed = str2double(get(controls.deltaspeed,'String'));
                delete(setwin);
            end
        end
        
%         function showerror(err)
%             warning(err.message);
%             for j = 1:size(err.stack,1)
%                 fprintf(2,'Error in ==> %s at %d\n',...
%                     err.stack(j).name,err.stack(j).line);
%             end
%         end
    end