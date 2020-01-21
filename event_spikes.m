function [spikes, spksevs, evspikes, foundindices, foundspikes, foundevspikes] =  event_spikes(allspikes, event, secbefore, secafter, parallel)
%evspikes =  event_spikes(allspikes, event, secbefore, secafter)
%returns subset of spike timestamps that occur
%  within time window around an event
%allspikes = timestamps of a given unit
%event = timestamps for one event
%secbefore, secafter = time (in secs) before and after each event
% OUTPUTS
% spikes: orig timestamps but only those around event
% spksevs: all the spike rates for each window
% evspikes: the spike ts as related to the event
% foundindices: the indices of the found spikes in a cell vector
% foundspikes: same, in a cell vector
%JRM
% JHB
% spkevs is the number of spikes that happened for each event over the
% event duration, much like sams spksevs

% added in foundindices to outputs, changed so that foundindices is saved
% for each event and its an array of cells for each event, shouldnt fuck up
% parent function.

% JHB 5-13-15
%allspikes=allspikes(:); event=event(:);
if ~exist('parallel','var')
    parallel=0;
end

numevents = length(event);

numfoundspikes = 0; % total number of spikes
spikes=[]; % vector of raw ts for each event
spksevs=[]; % vector of the # of spikes found for every event
foundindices={[]}; % cell mat of the indices per event
foundspikes={[]}; % cell mat of the raw spike ts per event
foundevspikes={[]}; % cell mat of found spikes locked to event
evspikes=[]; % vector of all the spikes locked to eventstarts

if length(secbefore)==1, secbefore=repmat(secbefore,length(event),1); end
if length(secafter)==1, secafter=repmat(secafter,length(event),1); end
if any(parallel)
    parfor ev = 1:numevents
        % find spikes and indices that are
        foundindices{ev} = find( (allspikes(:) > event(ev)-secbefore(ev)) & (allspikes(:) < event(ev)+secafter(ev)));
        foundspikes{ev} =allspikes( (allspikes(:) > event(ev)-secbefore(ev)) & (allspikes(:) < event(ev)+secafter(ev)));
        foundevspikes{ev}=foundspikes{ev}-event(ev); % lock to the trial start
        % take the number of spikes over your duration
        spksevs(ev)=length(foundindices{ev})/(secafter(ev)+secbefore(ev));
    end
else
    for ev = 1:numevents
        % find spikes and indices that are
        foundindices{ev} = find( (allspikes(:) > event(ev)-secbefore(ev)) & (allspikes(:) < event(ev)+secafter(ev)));
        foundspikes{ev} =allspikes( (allspikes(:) > event(ev)-secbefore(ev)) & (allspikes(:) < event(ev)+secafter(ev)));
        foundevspikes{ev}=foundspikes{ev}-event(ev); % lock to the trial start
        % take the number of spikes over your duration
        spksevs(ev)=length(foundindices{ev})/(secafter(ev)+secbefore(ev));
    end %end of for loop
end
spikes=cell2mat(foundspikes');
evspikes=cell2mat(foundevspikes');
numfoundspikes=sum(cellfun(@(a) length(a), foundindices'));
end
