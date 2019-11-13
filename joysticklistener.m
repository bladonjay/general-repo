classdef joysticklistener < handle
    % initialize with obj=joysticklistener()
    % two optional inputs: obj (random number) and period(default=10 ms)
    % Has properties:
    % id, joy, t, que, a, b, p, s
    
    % to query button presses
    % [events]=getEvents(joy)
    % this spits out a cell of variables first is button number second is
    % the time it got pressed
    
    properties (SetAccess = private, GetAccess = private)
        id
        joy
        t
        queue = zeros(0,3);
        a0
        b0
        p0
        s0
    end
    
    properties (Dependent)
        AveragePeriod
    end
    
    methods
        function obj = joysticklistener(id, period)
            if(nargin < 1 || isempty(id) || id == 0); id = 1; end
            if(nargin < 2 || isempty(period) || period == 0); period = 0.01; end
            
            obj.id = id;
            
            try obj.joy = vrjoystick(obj.id);
            catch ME
                error('BK:joysticklistener:notconnected',...
                    'Joystick with id ''%d'' not found.',obj.id);
            end
            
            [obj.a0, obj.b0, obj.p0] = read(obj.joy);
            obj.s0 = [obj.b0, obj.a0, obj.p0];

            obj.t = timer('Period',period,'ExecutionMode','fixedRate',...
                'Name','joysticklistenertimer',...
                'ErrorFcn',@joysticklistener.errorfcn);
        end
            
        function evs = getEvents(obj)
            evs = obj.queue;
            obj.queue = zeros(0,3);
        end
        
        function [a, b, p] = read(obj)
            [a, b, p] = read(obj.joy);
        end
        
        function check(obj, ts)
            [a, b, p] = read(obj.joy);
            deltaaup = (a~=obj.a0 & a <= -0.5);
            deltaadn = (a~=obj.a0 & a >= 0.5);
            deltab = (b~=obj.b0 & b ~= 0);
            deltap = (p~=obj.p0 & p ~= -1);

            delta = [deltab deltaaup deltaadn deltap];

            obj.a0 = a; obj.b0 = b; obj.p0 = p;
            obj.s0 = [b a a p];

            evs = find(delta)';
            evs(:,2) = ts;
            evs(:,3) = obj.s0(delta)';

            obj.queue = [obj.queue; evs];
        end
        
        function start(obj)
            if(isvalid(obj.t))
                obj.t.UserData = obj;
                obj.t.TimerFcn = ...
                    @(tobj,evt) check(tobj.UserData,datenum(evt.Data.time));
                start(obj.t);
            else error('BK:joysticklistener:invalidtimer',...
                    'Built in timer isinvalid.');
            end
        end
        
        function stop(obj)
            if(isvalid(obj.t))
                stop(obj.t);
                obj.t.TimerFcn = [];
                obj.t.UserData = [];
            end
        end
        
        function status = running(obj)
            if(isvalid(obj.t))
                status = strcmp(obj.t.Running,'on');
            else
                warning('BK:joysticklistener:invalidtimer',...
                    'Built in timer isinvalid.');
                status = false;
            end
        end
        
        function period = get.AveragePeriod(obj)
            period = obj.t.AveragePeriod;
        end
        
        function delete(obj)
            close(obj.joy);
            if(~isempty(obj.t) && isvalid(obj.t))
                stop(obj.t);
                delete(obj.t);
            end
        end
    end
    
    methods(Static)
        function errorfcn(t, evt)
            fprintf(2,'%s: %s\n',evt.Data.messageID,evt.Data.message);
            delete(t);
        end
    end
end
