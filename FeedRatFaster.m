function FeedRatFaster(a,PinMap,revs)
% function [ output_args ] = FeedRatFast(a,PinMap,revs,pacer)

% INPUTS:
%   a=          arduino object
%   PinMap=     Cell array of pin names, each pin name must be a string
%   revs=       Number of ticks to spin, is different for each motor
% OUTPUTS:
% feeds rat faster than the other one, on a stepper motor
%   Detailed explanation goes here


% pacer is to slow it down so other fx can execute need be



PinMap=[PinMap PinMap];
onoff=[1 0 1 0 0 1 0 1];
%onoff=[1 1 0 0 1 1 0 0];
regpos=repmat(1:8,1,revs);
postracker=1;

    
% each iteration advance
for i=1:length(regpos)
        writeDigitalPin(a,PinMap{regpos(postracker)},onoff(regpos(postracker)));
        postracker=postracker+1;
        %fprintf('test \n');
end

% now turn all off
for i=1:length(PinMap)
    writeDigitalPin(a,PinMap{i},0);
end



end