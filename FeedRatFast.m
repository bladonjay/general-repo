function [mytimer] = FeedRatFast(a,PinMap,revs,pacer)
% function [ output_args ] = FeedRatFast(a,PinMap,revs,pacer)

% INPUTS:
%   a=          arduino object
%   PinMap=     Cell array of pin names, each pin name must be a string
%   revs=       Number of ticks to spin, is different for each motor
%   pacer=      This slows down the spin rate if you want to use other fx
% OUTPUTS:
%   mytimer=    This is the Timer object
%feeds rat faster than the other one, on a stepper motor
%   Detailed explanation goes here


% pacer is to slow it down so other fx can execute need be
if ~exist('pacer','var')
    pacer=.02;
end


mytimer=timer('name','FeederTimer','TimerFcn',@timerfcn,...
    'StopFcn',@stopfcn,'Period',pacer,'ExecutionMode','FixedDelay',...
    'TasksToExecute',revs*8,'ErrorFcn',@Errfcn,'BusyMode','queue');

PinMap=[PinMap PinMap];
onoff=[1 0 1 0 0 1 0 1];
%onoff=[1 1 0 0 1 1 0 0];
regpos=repmat(1:8,1,revs);
postracker=1;
start(mytimer);
    
    % each iteration advance
    function timerfcn(varargin)
        writeDigitalPin(a,PinMap{regpos(postracker)},onoff(regpos(postracker)));
        postracker=postracker+1;
        %fprintf('test \n');
        % now kill all pins so motor doesnt heat up
    end


    function stopfcn(varargin)
        for i=1:length(PinMap)
            writeDigitalPin(a,PinMap{i},0);
        end
    end

    function Errfcn(varargin)
        for i=1:length(PinMap)
            writeDigitalPin(a,PinMap{i},0);
        end
        delete(mytimer);
    end

end