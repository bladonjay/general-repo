function [spike_times] = pe_raster_endsample_phase(spikes, events, secbefore, secafter, end_interval, short_pre_sample,phase)
%function pe_raster(spikes, events, secbefore, secafter,color)
% plots a perievent raster
% spikes=array of timestamps for a single unit
% events=array of timestamps for a single behavioral event
% secbefore=time (in seconds) before event
% secafter=time (in seconds) after event
% color can be a string (e.g., 'b') or a rgb triplet (e.g., [0 0 1])
% spikes_times = spikes times relative to events
% **calls event_spikes.m
% phases is a companion matrix of the phase of each spike in the LFP
%JRM
num_events = length(events);
hold on
spike_times = [];
for ev_num = 1:num_events
    %find the spikes within a time window around each event
    [tempspikes,~,phaseindices] = event_spikes(spikes, events(ev_num),secbefore, secafter);
    tempphases=phase(phaseindices{1});
    % have to fetch those phases too
    
    %find the offset (in secs) from the current event
    found = 0;
    if (tempspikes > 0)
        found = 1;
        tempspikes = tempspikes - events(ev_num);
        %temp_end_interval = end_interval(ev_num) - events(ev_num);
        spike_times = [spike_times tempspikes];
        lineheight = [zeros(1,length(tempspikes)); .95*ones(1,length(tempspikes))]      + ev_num;
        
        if ~isempty(phase)
            % phases come in from -pi to pi so:
            tempphases=abs(tempphases);
            colormap_fr =linspace(0,pi,255);
            [~,colormap]=histc(tempphases,colormap_fr);
            rgbcolormap=ind2rgb(colormap,jet(256)); % should use hsv so its circular
            % now do each line with its color phase
            for q=1:length(tempspikes)
                %if tempphases(q)<pi/2
                linehandle = line([tempspikes(q); tempspikes(q)], lineheight);
                set(linehandle,'LineWidth',2, 'color', rgbcolormap(q,1,:));
                % just cause i like to see my progress
                %drawnow
                hold on
                %end
            end
        else
            linehandle = line([tempspikes; tempspikes], lineheight);
            set(linehandle,'LineWidth',2, 'color', 'k');
        end

% Black sample end marks removed temporarily CK 5/16/13        
%         temp_end_interval = end_interval(ev_num) - events(ev_num);
%         lineheight = [zeros(1,length(temp_end_interval)); .95*ones(1,length(temp_end_interval))]      + ev_num;
%         linehandle = line([temp_end_interval; temp_end_interval], lineheight);
%         set(linehandle,'LineWidth',2, 'color', 'k');
%         
%         if (short_pre_sample(ev_num) > 0)
%             temp_pre_sample = short_pre_sample(ev_num) - events(ev_num);
%             linehandle = line([temp_pre_sample;  temp_pre_sample], lineheight);
%             set(linehandle,'LineWidth',2, 'color', 'k');
%         end
%     else
%         if (short_pre_sample(ev_num) > 0)
%             temp_pre_sample = short_pre_sample(ev_num) - events(ev_num);
%             lineheight = [zeros(1,length(temp_pre_sample)); .95*ones(1,length(temp_pre_sample))]      + ev_num;
%             linehandle = line([temp_pre_sample;  temp_pre_sample], lineheight);
%             set(linehandle,'LineWidth',2, 'color', 'k');
%         end
%         temp_end_interval = end_interval(ev_num) - events(ev_num);
%         lineheight = [zeros(1,length(temp_end_interval)); .95*ones(1,length(temp_end_interval))]      + ev_num;
%         linehandle = line([temp_end_interval; temp_end_interval], lineheight);
%         set(linehandle,'LineWidth',2, 'color', 'k');
%         
    end%end of if found
end%end of for loop

xlim([secbefore*-1 secafter]);
ylim([0 num_events+1]);

hold off

