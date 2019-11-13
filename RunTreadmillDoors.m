function RunTreadmillDoors(settings)

% $Id: RunTreadmill.m 4831 2013-05-13 19:56:09Z nrobinson
% commented by JH Bladon
% RunTreadmill runs a treadmill task with a gui, will open a treadmill
% object, a joystick object (if possible) and will try to connect to the
% plexon map server.

% function showboxes and function mazeinfo_treadmill are where you change
% your boxes for triggering events automatically





% set time stream
rstream = RandStream('mt19937ar','Seed', mod(prod(clock()),2^32));
% set the global time stream
RandStream.setGlobalStream(rstream);


        
defsettings.laserpin=2; % lazer pin is 2 (nicks choice)
defsettings.usedoor=true; % default is to use the door
defsettings.lightpin=4; % which pin is for the delay light
defsettings.door1pin=3;  % default door pin is 3
defsettings.door2pin=5;  % default door pin is 3
defsettings.door1bell=7;
defsettings.door2bell=8;
defsettings.tmbell=12; % because 9 10 and 11 are reserved for pwm
defsettings.tmsbell=13;
defsettings.door1down=.55; %how far down to put the door
defsettings.door2down=.55; %how far down to put the door

defsettings.doorwait=2; % when to put the door back down



% a bunch of default settings
defsettings.mapserver = false;       % to sync up to map server?

defsettings.joystick = true;         % use a joystick
defsettings.tracking = false;        % dont use tracking a trigger
defsettings.maze = 'TreadmillJ';      % its a treadmill maze
defsettings.training= false;         % not a training session
defsettings.deltaspeed = 5;          % change speed?
defsettings.laserdelay = 2;          % seconds after tmill start for laser on
defsettings.laserduration = 2;       % seconds
defsettings.rundelay = 0;            % seconds after door that treadmill will run
defsettings.quickdelay=5;            % this is for the indicator light for sampling
defsettings.lasermode = 'nolaser';   % 'nolaser' = no laser
                                     % '5050'    = 50% no laser, 50% laser
                                     % 'arm'     = 20% arm, 40% laser, 40% no laser
 
defsettings.treadmill.on = true;     % turn treadmill on?
defsettings.treadmill.speeds = 30;   % at 30 cm/sec
defsettings.treadmill.times = 2 ;   % with 2 sec delay
defsettings.treadmill.distances = 1500;  % we dont need either
defsettings.treadmill.fixed = 'time'; % Fixed parameter
defsettings.treadmill.vary = 'speed'; % Independant parameter

defsettings.treadmill.pausenotstop = true; % Remember settings when manually stopping treadmill


%% initialize treadmill, timer, and arduino
% initiate treadmill object in the com 4 port
tr= Treadmill('COM4');
% initiate serial object for lazer and doors (arduino uno)
%laserserial = serial('COM5','BaudRate',9600);
% instead of using serial, use the laser object
myarduino=arduino();
door1=servo(myarduino,defsettings.door1pin);
door2=servo(myarduino,defsettings.door2pin);


% set up our time recorder, that will use a 0.01 second update rate
tim = timer('StartFcn',@startfcn,'TimerFcn',@timerfcn,...
    'StopFcn',@stopfcn,'ErrorFcn',@errorfcn,'Period',0.01,...
    'ExecutionMode','fixedDelay','Name','Main RunTreadmill Timer');

% plex server address(start at zero, it will find the address on startup)
plexserver = 0;

% dont use device num, but we could
deviceNum = 0;

% do you need to save on shutdown?
needtosave = false;

EVENTCODES = maze_events();       % Load the event codes. gotta add a door code
% start a history seed
history = [];

% start the timer
timeref = tic;


%fopen(laserserial);

% make sure you close all objects
%autocloselaserserial = onCleanup(@() fclose(laserserial));
autocloselaserserial = onCleanup(@() delete(instrfindall));

%% session annotations
lap=1; % start at lap 1
sessionnum = 0;             % not crucial but keeps track of the session number
session.comment = '';       % no default comments yet
session.date = now();       % all annotations
session.laps = [];          % total laps
session.correct = [];       % also something to add
session.dig = [];           %
session.objects = [];       %
session.laser = [];         %
session.samples = [];       %

trialblocksize = 20;        % this will designate block size pre randomization
randomizer = [];
objectlist = [];
laserlist = [];

% these are for other tasks
treadmillstart = 0;
treadmillcount = 0;
treadmillind = 0;
treadmillcurspeed = 0;
treadmillcurtime = 0;
treadmillpausestate = 0;

% start with tracking off/not triggered
trackstat.ready = false;
% joy is joystick (off for now)
joy = 0;

%% initiate everything with default settings
if(nargin == 0)
    settings = checksettings();
else
    settings = checksettings(settings);
end
savetrackstat = trackstat; % DEBUGGING


%% open the gui
h = waterportsGUI();
% Set close controls
set(h.controls,'CloseRequestFcn',@closeButton)
set(h.timer,'CloseRequestFcn',@toggleTimer);

% now set callback functions
% set(h.c,'Callback',{@(x,y) process_event(EVENTCODES.,'computer')});
% set(h.i,'Callback',{@(x,y) process_event(EVENTCODES.,'computer')});
set(h.udoor,'Callback',{@ movedoor,door2,settings.door2down,settings.door2bell} ); % default is to use the door 
set(h.mdoor,'Callback',{@ movedoor,door1,settings.door1down,settings.door1bell});
set(h.tmillstart,'Callback',{@(x,y) process_event(EVENTCODES.TM,'computer')});
set(h.tmillstop,'Callback',{@(x,y) process_event(EVENTCODES.TMS,'computer')});
set(h.open,'Callback',@openHist);
set(h.save,'Callback',@saveHist);
set(h.pause,'Callback',@playpause);
set(h.reset,'Callback',@resetchecker);
set(h.settings,'Callback',@settingscallback);
set(h.showtimer,'Callback',@toggleTimer);
applysettings;


%% timer functions
    function startfcn(varargin)
        output(h.maindsp, 'Disabling Screen Saver...');
        [status, result] = dos('REG ADD "HKCU\Control Panel\Desktop" /v ScreenSaveActive /t REG_SZ /d 0 /f');
        if(status == 0); output(h.maindsp, 'Screen Saver Disabled.');
        else output(h.maindsp, regexprep(result,'\n',' '));
        end
        
        % Open connection with the Plexon Server
        
        % variable plexserver is the server ID used for al plx functions
        if(settings.mapserver)
            output(h.maindsp,'Connecting to Plexon MAP Server...');
            try plexserver = PL_InitClient(0);
            catch %#ok<CTCH>
                plexserver = 0;
            end
            
            if plexserver == 0
                output(h.maindsp,'Failed to connect to Plexon MAP Server');
            end
        end
        
        if(plexserver ~= 0)
            output(h.maindsp, 'Connected to Plexon MAP Server.');
            
            pars = PL_GetPars(plexserver);
            if(~pars(11))
                output(h.maindsp, 'Sort Client is not running.');
                warndlg('Sort Client must be running to use the button box or video tracking. Disconnecting from MAP Server.','Sort Client not running.');
                output(h.maindsp, 'Closing connection to Plexon MAP Server...');
                PL_Close(plexserver);
                plexserver = 0;
            end
        end
        
        set(h.pause,'String','Pause');
        set(h.pause,'Enable','on');
        

    end

    
    function errorfcn(varargin)
        output(h.maindsp, 'Error encountered.');
    end

    % stops the timer, which is really the pase
    function stopfcn(varargin)
        output(h.maindsp, sprintf('Pausing...'));
        
        if (plexserver ~= 0);
            output(h.maindsp, 'Closing connection to Plexon MAP Server...');
            PL_Close(plexserver);
            plexserver = 0;
        end
        
        output(h.maindsp, 'Enabling Screen Saver...');
        [status, result] = dos('REG ADD "HKCU\Control Panel\Desktop" /v ScreenSaveActive /t REG_SZ /d 1 /f');
        if(status == 0); output(h.maindsp, 'Screen Saver Enabled.');
        else output(h.maindsp, regexprep(result,'\n',' '));
        end
        
        set(h.pause,'String','Start');
        set(h.pause,'Enable','on');
    end



% timer function will execute every 0.01 seconds once you press the
% start button on the GUI
    function timerfcn(varargin)
        % if we have a joystick(settings), lets see if it says anything
        % what number is joy
        if(settings.joystick && joy ~= 0 && isvalid(joy))
            joyevents = getEvents(joy);
            for k = 1:size(joyevents,1)
                % for however many buttons are pressed, process that event
                joystick_event(joyevents(k,1),joyevents(k,2));
            end
        end
        % if were using plexon,
        if(plexserver ~= 0)
            % not sure why we have to wait ten seconds
            res = PL_WaitForServer(plexserver, 10);
            if(res ~= 0)
                [~, t] = PL_GetTS(plexserver);
                
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


%% % joystick and events
    function joystick_event(event_num,event_time)
        % bladon map for controller
        jsmap=SNESgamepadmap;
        %jsmap = joystickmap;
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
            case EVENTCODES.TM % treadmill start (101)
                output(h.maindsp, sprintf('Starting treadmill (%.2f)', sessionts));
                initTreadmill(event_num,event_time);
                if(plexserver ~= 0); PL_SendUserEvent(plexserver, event_num); end
                
            case EVENTCODES.TMS % treadmill stop (108)
                output(h.maindsp, sprintf('Stopping treadmill (%.2f)', sessionts));
                stopTreadmill();
                if(plexserver ~= 0); PL_SendUserEvent(plexserver, event_num); end
                
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
                if(plexserver ~= 0); PL_SendUserEvent(plexserver, event_num); end
                
            case {EVENTCODES.DIG, EVENTCODES.NODIG} % Responses (105 or 106)
                ratresponse(event_num, event_time)
                if(plexserver ~= 0); PL_SendUserEvent(plexserver, event_num); end
                
            case EVENTCODES.LASERZONE % 101, 
                armlaser(event_num, event_time)
                if(plexserver ~= 0); PL_SendUserEvent(plexserver, event_num); end
                
                
            case EVENTCODES.DOOR1
                movedoor(1,1,door1,settings.door1down,settings.door1bell)
                if(plexserver ~= 0); PL_SendUserEvent(plexserver, event_num); end
                
            case EVENTCODES.DOOR2
                movedoor(1,1,door2,settings.door2down,settings.door2bell)
                if(plexserver ~= 0); PL_SendUserEvent(plexserver, event_num); end
            
            case EVENTCODES.QuickTimer
                qt=timer('StartFcn',@lighton,'TimerFcn',@lightoff,'StartDelay',...
                    settings.quickdelay,'TasksToExecute',1);
                %temp= @(pin) writeDigitalPin(a,pin,1);
                start(qt);
                
        end
        updateCounter();
        updateDisplay();
    end



%% tabulating the trial types and rat responses
% i'll need to revamp this so that I can add in indicator buttons
    function createtrialtypes()
        if(settings.laserduration <= 0 && strcmpi(settings.lasermode,'5050'))
            settings.lasermode = 'nolaser';
        end
        
        switch settings.lasermode
            case 'arm' % return arm lazer
                % permute 1:20 and pull evens (thats obj list)
                randomizer = mod(randperm(trialblocksize),20);
                objectlist = (mod(randomizer,2)==0); % take even trials
                
                laserlist = nan(size(randomizer));
                laserlist(randomizer<11) = 0;       %0 % no laser change back for old 'arm'
                laserlist(randomizer>=11) = 2;      %2 % treadmill laser (second half)
                % laserlist(trialtypes>=16) = 2;     % arm laser (16:20)
            case '5050' % half trials are lazer
                randomizer = mod(randperm(trialblocksize),20);
                objectlist = (mod(randomizer,2)==0);
                laserlist = (randomizer>9); % laserlist is trials above 9
            case 'allLaser'
                randomizer = mod(randperm(trialblocksize),20);
                objectlist = (mod(randomizer,2)==0);
                laserlist = ones(size(randomizer));
            otherwise % 
                randomizer = mod(randperm(trialblocksize),20);
                % for no objects, only correct or incorrect response
                % zero means dig
                objectlist=zeros(size(randomizer));
                % for two objects
                %objectlist = (mod(trialtypes,2)==0);
                laserlist = false(size(randomizer));
        end
        
        % this makes sure there are no more than 3 repeats in a row
 %        [objectlist]=trialGenCM(objectlist');
    end


%
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

%% lazer, door and treadmill functions
    
    % arm lazer is basically the wrapper to make sure the lazer should be
    % going, it arms the lazer and checks to see if its the correct trial
    % type for the lazer to fire
    function armlaser(event_num, event_time)
        lap_ = mod(lap-1,numel(objectlist))+1;
        if(~isempty(laserlist) && lap > 0 && laserlist(lap_)==2 &&...
                settings.laserduration > 0)
            %pulselaser(laserserial,0,settings.laserduration);
            
            % pulse the lazer, on for x seconds...
            pulselaser2(myarduino,settings.laserpin,0,settings.laserduration);
            
            output(h.maindsp, 'Pulsing Laser (arm)');
        end
    end

% initTreadmill starts the treadmill sequence, which could include a
% waterport opening, followed by a pause, followed by the actual start
% of the treadmill.
% in this case its just start treadmill
    function initTreadmill(event_num, event_time)
        if(strcmp(tr.Running,'on'))
            output(h.maindsp, 'Treadmill is already running.');
        else
            needtosave = true;
            tr.SetSpeed=settings.treadmill.speeds;
            tr.start(settings.treadmill.times) %+settings.rundelay)
            
            writeDigitalPin(myarduino,settings.tmbell,1); writeDigitalPin(myarduino,settings.tmbell,0); 
            
            door1timer=pulseservo(myarduino,door1,settings.treadmill.times,settings.doorwait,settings.door1down,settings.door1bell);
            lap_ = mod(lap-1,numel(objectlist))+1;
            
            if(~isempty(laserlist) && lap > 0 && laserlist(lap_)==1 &&...
                    settings.laserduration > 0)
                pulselaser2(myarduino,settings.laserpin,settings.laserdelay,settings.laserduration);
                %pulsearduino(myarduino,settings.laserpin,settings.laserdelay,settings.laserduration);
                
                
                output(h.maindsp, 'Pulsing Laser (treadmill)');
            end
        end
    end

% stopTreadmill is called to force the treadmill to stop immediately.
    function stopTreadmill()
        % make sure door wont open
       
        
        % stop treadmill
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
        stopdooring();
        writeDigitalPin(myarduino,settings.tmsbell,1); writeDigitalPin(myarduino,settings.tmsbell,0);
    end

% change door changes the position of the door
    
    function movedoor(~,~,door,down,bell)
        stopdooring();
        isup=readPosition(door);
        if isup==0
            writePosition(door,down)
        else
            writePosition(door,0);
        end
        
        writeDigitalPin(myarduino,bell,1); writeDigitalPin(myarduino,bell,0); 
        
    end

    % dont use
%     function laseroff(varargin) 
%         writeDigitalPin(myarduino,settings.laserpin,0);
%     end
    function stopdooring()
        stopdoor=timerfind('name','doortimer');
        deets=whos('stopdoor');
        if strcmpi(deets.class,'timer')
            output(h.maindsp,'found door timer, stopping door')
            stop(stopdoor);
            %delete(stopdoor);
        
        else
            %output(h.maindsp,'cant get door timer');
        end    
    end

 % quick light indicator functions
    function lighton(varargin)
        thists = (now() - session.date)*24*60*60;
        writeDigitalPin(myarduino,settings.lightpin,1);
        output(h.maindsp, sprintf('activated quick timer(%.2f)',thists));
    end

    function lightoff(varargin)
        thists = (now() - session.date)*24*60*60;
        writeDigitalPin(myarduino,settings.lightpin,0);
        output(h.maindsp, sprintf('deactivated quick timer(%.2f)',thists));
        
    end
%% for modifying settings (opens new window)
    function settingscallback(varargin) %#ok<VANUS>
        settings = changesettings(settings, h);
        settings = checksettings(settings);
        applysettings;
    end

    function applysettings()

        
        % set treadmill settings, if this maze has a treadmill, enable the
        % treadmill functions
        if strfind(regexprep(lower(settings.maze),'[^a-z]',''),'treadmill')
            set(h.tmillstart,'Enable','on');
            set(h.tmillstop,'Enable','on');
        else
            set(h.tmillstart,'Enable','off');
            set(h.tmillstop,'Enable','off');
        end
        
        if(isempty(laserlist) || ...
                (any(laserlist==2) && ~strcmpi(settings.lasermode,'objectlaser')) || ...
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
            if(joy == 0 || ~isvalid(joy));
                try
                    joy = joysticklistener(1);
                    start(joy);
                catch %#ok<CTCH>
                    joy = 0;
                    warndlg('Joystick failed to initialize.','Joystick Disabled');
                    settings.joystick = false;
                    
                end
            end
        elseif(joy ~= 0)
            if(isvalid(joy)); delete(joy); end
            joy = 0;
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
% check fieldtypes makes sure that new struct doesnt have missing
% or messed up fields (it only turns bad fields into old struct
% fields, othewise it takes the new structs fields
    function newstruct = checkfieldtypes(newstruct, refstruct)
        fields = fieldnames(refstruct);
        
        for j = 1:size(fields,1)
            % if new struct doesnt ahave the field, or its class is
            % wrong, new struct takes old structs field
            if(~isfield(newstruct,fields{j}) || ...
                    ~isa(newstruct.(fields{j}),class(refstruct.(fields{j}))))
                newstruct.(fields{j}) = refstruct.(fields{j});
                % if newstruct has a field but its empty, or it has a nan
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
%% saving and history functions
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
%% Update and display functions
    function updateDisplay()
        % add first name into settingsstring, like laser: major, and
        % doors: guillatine
        settingsstr = sprintf('Doors');
        
        % if we have treadmill on put it in the display
        if(settings.treadmill.on); settingsstr = [settingsstr ', Treadmill: on'];
        else settingsstr = [settingsstr ', Treadmill: off'];
        end
        % update the display
        set(h.settingsdsp,'String',settingsstr);
        
        % now update the lap info
        if(lap < 0)
            sessionstr = sprintf('Not Started');
        elseif(lap == 0)
            sessionstr = sprintf('Lap: %d, T1: %g, T2: %g',lap,tim.AveragePeriod,joy.AveragePeriod);
        else
            sessionstr = 'Lap: %d, Correct: %d (%.2f%%), T1: %g, T2: %g';
            sessionstr = sprintf(sessionstr,lap,sum(session.correct),...
                sum(session.correct)/lap*100,tim.AveragePeriod,joy.AveragePeriod);
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

% this is for the main display, not sure why he has to do this...
    function output(display, str)
        output_list = get(display,'String');
        output_list = output_list(end:-1:1);
        num = length(output_list)+1;
        str = strtrim(str);
        output_list{end+1} = sprintf('%3d  %s',num,str);
        set(display,'String',output_list(end:-1:1));
        drawnow
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
                % if youre using laser on the return arm
            elseif(laserlist(lap_)==2); lstr = 'L-Arm';
            elseif(laserlist(lap_)==3); lstr = 'L-Obj';
            else lstr = ' ';
            end
            
            % Display lap number, object, laser, and some other string...
            % residual from ben...
            set(h.counter,'String',sprintf('lap %3d: %1d %1s %2s',lap,objectlist(lap_),lstr));
        end
    end

    function updatelights()
        % not real function yet, but will be an indicator light for trial, and laser indicator
    end
%% reset the error checker
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
        
        randomizer = [];
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
%%
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
            if(joy ~= 0 && isvalid(joy)); delete(joy); end
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
            delete(instrfind);
            
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
                %maximize(h.timer);
        end
    end
%% Organize the windows for the gui thats made at top
    function h = waterportsGUI()
        screen = get(0,'ScreenSize');
        winsize = [450 600];
        
        h.controls = figure('Name','Treadmill Task',...
            'NumberTitle','off','MenuBar','none',...
            'Position', [10 (screen(4)-winsize(2)-30) winsize]);
        
        h.buttons = uipanel('Parent',h.controls,...
            'Units','normalized','Position',[0.01 0.82 0.35 0.17]);
        
        % looks like this is for the correct and incorrect choices
        %             h.c = uicontrol(h.buttons,'Style','pushbutton',...
        %                 'Units','normalized','Position',[0.02 0.35 0.47 0.30],...
        %                 'String','Correct');
        %             h.i = uicontrol(h.buttons,'Style','pushbutton',...
        %                 'Units','normalized','Position',,...
        %                 'String','Incorrect');
        
        % replacing the choice buttons with a door up/down button
        h.mdoor = uicontrol(h.buttons,'Style','pushbutton',...
            'Units','normalized','Position',[0.02 0.35 0.47 0.30],...
            'String','Move Door 1');
        
        
        h.udoor = uicontrol(h.buttons,'Style','pushbutton',...
            'Units','normalized','Position',[0.51 0.35 0.47 0.30],...
            'String','Move Door 2');
        
        
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
            'Position', [5 39 screen(3)-1000 screen(4)-800]);
        % e.g. the timer visual
        h.counter = uicontrol(h.timer,'Style','text',...
            'Units','normalized','Position',[0 0.0 1 .5],...
            'BackgroundColor','black','ForegroundColor','white',...
            'FontUnits','normalized','FontSize',.65);
        
        h.timerstr = uicontrol(h.timer,'Style','text',...
            'Units','normalized','Position',[0 0.5 1 .5],...
            'BackgroundColor','black','ForegroundColor','white',...
            'FontUnits','normalized','FontSize',.95);
    end

%% change settings window
    function settings = changesettings(settings, h)
        % nick has code to add in object laser
        LASERMODES = {'nolaser','5050','arm','allLaser'};
        
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
            'BackgroundColor','white','String',num2str(settings.laserduration));
        
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
        
        labels.rundelay = uicontrol(setwin,'Style','text',...
            'String','Run Delay:');
        
        controls.rundelay = uicontrol(setwin,'Style','edit',...
            'BackgroundColor','white','String',num2str(settings.rundelay));
        
        labels.timerdelay=uicontrol(setwin,'Style','text',...
            'String','Sample Timer');
         
        controls.timerdelay = uicontrol(setwin,'Style','edit',...
            'BackgroundColor','white','String',num2str(settings.quickdelay));       
        
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
            settings.rundelay = str2double(get(controls.rundelay,'String'));
            settings.timerdelay=str2double(get(controls.timerdelay,'String'));
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