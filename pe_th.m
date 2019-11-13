function [xtimes,rate_per_bin,maxrate] = pe_th(spikes, events, secbefore, secafter, nbins)
%function pe_th(spikes, events, secbefore, secafter,nbins)
%plots a perievent time histogram
%spikes=array of timestamps for a single unit
%events=array of timestamps for a single behavioral event
%secbefore=time (in seconds) before event (positive number)
%secafter=time (in seconds) after event (positive number)
% calling it looks like peth(units(1).ts,session.bf.XRE,2,2,10,5)
%**calls event_spikes.m
%JRM 5-27-05


%%
tot_duration = secbefore + secafter;
bin_duration = tot_duration/nbins;

num_events = length(events);
allfoundspikes = [];



for ev_num = 1:num_events
	%find the spikes within a time window around each event (centered
	%around that event)
	[tempspikes,~,tryspikes] = event_spikes(spikes, events(ev_num),secbefore, secafter);
	if any(tempspikes)
		%add those spikes to our accumulating total
		allfoundspikes = [allfoundspikes; tryspikes];
	end%end of if found
	
end%end of for loop

%determine the midpoint values for the histogram bins
for b = 1:nbins
	bin_times(b) = (-1*secbefore) +(.5*bin_duration) + (b-1)*bin_duration;
end

%get the counts per bin, but don't plot yet
[spikes_per_bin, xtimes] = hist(allfoundspikes,bin_times);
%convert spikes per bin to hertz
rate_per_bin = (spikes_per_bin ./ bin_duration) / num_events;
maxrate=max(rate_per_bin);
%plot it
%bar(xtimes,rate_per_bin,1);
%xlim([secbefore*-1 secafter]);%make sure the x scale is what we expect
%ylim([0 20]) %Dont care about rate cap because we will standardize
% later


end
