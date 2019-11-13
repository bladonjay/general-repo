function [mytimer] = FeedRat2(a,PinMap, revs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% posregister= [0 0 1 1 1 0 0 0;...
%               1 0 0 0 0 0 1 1;...
%               1 1 1 0 0 0 0 0;...
%               0 0 0 0 1 1 1 0];
%           
 posregister= [0 1 1 0 0 1 1 0;...
               1 0 0 1 1 0 0 1;...
               1 1 0 0 1 1 0 0;...
               0 0 1 1 0 0 1 1];
                   

mytimer=timer('name','BeamBreakTimer','TimerFcn',@timerfcn,...
    'StopFcn',@stopfcn,'Period',.05,'ExecutionMode','SingleShot',...
    'BusyMode','queue');

start(mytimer);

    function timerfcn(varargin)
            regpos=repmat(1:8,1,revs);


        for i=1:length(regpos)
            % to get multiple of 8
            %if regpos==0, regpos=8; end
            
            % cycle throguh to move gears
            writeDigitalPin(a,PinMap{1},posregister(1,regpos(i)));
            writeDigitalPin(a,PinMap{2},posregister(2,regpos(i)));
            writeDigitalPin(a,PinMap{3},posregister(3,regpos(i)));
            writeDigitalPin(a,PinMap{4},posregister(4,regpos(i)));
        end
    end

    function stopfcn(varargin)
        delete(mytimer);
    end
end

