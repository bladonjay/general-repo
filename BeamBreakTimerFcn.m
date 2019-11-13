function [mytimer,sensorvals]=BeamBreakTimerFcn(arduino,pinreg,revs)

mytimer=timer('name','BeamBreakTimer','TimerFcn',@timerfcn,...
    'StopFcn',@stopfcn,'Period',.03,'ExecutionMode','FixedSpacing',...
    'BusyMode','queue','StartDelay',.2);




%figure;
sensorv=1;
sensorvals=[];
start(mytimer);





    function timerfcn(varargin)

            regpos=repmat(1:8,1,revs);
            onoff=[1 0 1 0 1 0 1 0];
            for i=1:50
                writeDigitalPin(arduino,pinreg{regpos},onoff(regpos));
            end
            
            stop(mytimer)
    end
    function stopfcn(varargin)
        delete(mytimer);
        
    end
end
