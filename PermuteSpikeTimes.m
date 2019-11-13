function [newspiketimes,bracketedspiketimes] = PermuteSpikeTimes(spiketimes, coords, minshift)
% Function PermuteSpikeTimes randomly shuffles spiketimes while keeping
% the temporal order. coords are to make sure you have the session length


% it basically performs a circular shift starting at a random timepoint
% there is a border window so that it doesnt shift all spikes by 1 though

% first get spike times that are within our coordinates (becaues we have
% already epoched the timestamps)
% log our old ones
newspiketimes=[];
% in seconds
if ~exist('minshift','var')
    minshift=20;
elseif isempty(minshift)
    minshift=20;
end


% find start and end of epoch
start=min(coords(:,1)); finish=max(coords(:,1));

% take only spikes within that interval
spiketimes2=spiketimes(spiketimes>=start & spiketimes<finish);

% transform spiketimes to start at 0
spiketimes3=spiketimes2-start;

% make new strobe data: continuous at 100 hz with border (minshift)
% this cuts ends off too for a danger zone to not shift too little
cut_pool=minshift: .1 :finish-start-minshift;

% if we have enough spikes for analysis
%if length(spiketimes3)>2
    
    % now pull random timestamp from numbergenerator
    [cutposition]=datasample(cut_pool,1);
    

    
    % now circular shift those numbers
    
    % all numbers below our point get shifted forward (Basically add our
    % newstart to these data
    firsthalf=spiketimes3(spiketimes3<cutposition);
    
    % add the end chunk on front
    newsecondhalf=firsthalf+(finish-cutposition);
    
    % all numbers above our point start from 0 (secondhalf first ind will be 0)
    secondhalf=spiketimes3(spiketimes3>=cutposition);
    
    % now subtract our cutposition from those, and then put them together
    % with the new second half
    shuffledspikes=[secondhalf-cutposition; newsecondhalf];
    
    % now put the spikes back to the start time
    newspiketimes=shuffledspikes+start;
    
    % and save our bracketed spike times (original times in that epoch)
    bracketedspiketimes=spiketimes2;
%else
    %newspiketimes=nan;
end

