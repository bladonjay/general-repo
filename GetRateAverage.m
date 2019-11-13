function [averagerate,session_duration] = GetRateAverage(SpikeTimes,trackingts,jumptolerance)
% function [average rate, session duration]= GetRateAverage(SpikeTimes,
% Tracking timestamps, jumptolerance
%Gets the session average firing rate by using only epochs where tracking
%was continuous, or, if you epoch your data it will take the rate only
%during your epochs

% INPUTS
%   SpikeTimes = the spike times of your cell (vector of timestamps)
%   Tracking  = an epoch of tracking data (also a vector of timestamps)
%   Jumptolerance = smaallest gap between tracking ts that youre ok combing
%   over, this is important for trial epoched data
% OUTPUTS
% Average rate is the average firing rate across your epoch
% Session Duration is the duration  of time that the average was calculated
% over, also a surrogate for the total time spent in that 'epoch'


if ~exist('jumptolerance','var');
    jumptolerance=1;
elseif isempty(jumptolerance)
    jumptolerance=1;
end



% JH Bladon, 5-27-15

% first have to make sure dimensions are all good:
temp=SpikeTimes; clear SpikeTimes;
SpikeTimes(:,1)=temp; clear temp;

temp=trackingts; clear trackingts;
trackingts(:,1)=temp; clear temp;


dt = diff(trackingts); %difference between consecutive coords

% not sure why i need to redo coords
%$ tcoord = trackingts((2:end),1); % remove coords so dt1 is for before tcoord1

% again Im going to just use the original tracking ts here
% betweens is the midpoint between each tracking timestamp
%$  betweens=trackingts(2:end,1)-(dt/2); % now betweens is just before tcoords
% this is the time window that continuous time data will be assigned that
% tcoord

% %betweens(end+1)=trackingts(end,1);
% first remove gaps in tracking

timejumps = [];
% timejump values are the start of the gap

timejumps = cat(1,timejumps,find(diff(trackingts) > jumptolerance));
% data isnt interpolated, so we will have to use a timejump window of
% around 10, not .334 (which is strobe rate)

% each timejump is the interval of that timestamp to the one AFTER it
EpochedSpikes = [];
session_duration = 0;
% log session duration and pull spikes that are ONLY within our epochs,
% epochs, that is, that have ocntinuous tracking (diff is less than 10
% second


% this only remove spikes, this doesnt modify any coordinates.  It also
% logs the duration of each epoch i keep and sums them at the end
if length(timejumps) >= 1
    for i = 1:length(timejumps)+1
        if i == 1 % if its the first epoch
            % remove ts that will be assigned earlier than ts1 and after
            % the ts that starts your first jump
            EpochedSpikes = cat(1,EpochedSpikes,SpikeTimes(SpikeTimes >= trackingts(1) & SpikeTimes <= trackingts(timejumps(i))));
            session_duration = trackingts(timejumps(i)) - trackingts(1);
        elseif  i == length(timejumps)+1 % if its the last epoch
            % remove the ts that will be assigned to the first ts of the
            % last epoch
            % start one timestamp after the last pause (or the end of the
            % gap)
            EpochedSpikes = cat(1,EpochedSpikes,SpikeTimes(SpikeTimes >= trackingts((timejumps(i-1)+1)) & SpikeTimes <= trackingts(end)));
            session_duration = session_duration + trackingts(end) - trackingts(timejumps(i-1)+1);
        else % middle epochs
            % remove spikes that will be assigned to the ts that border
            % the pause (each will be assigned a large lapsed time
            EpochedSpikes = cat(1,EpochedSpikes,SpikeTimes(SpikeTimes >= trackingts((timejumps(i-1)+1)) & SpikeTimes <= trackingts(timejumps(i))));
            session_duration = session_duration + trackingts(timejumps(i)) - trackingts(timejumps(i-1)+1);
        end
    end
    % would be nice to get these in ascending order huh
    SpikeTimes = sort(EpochedSpikes);
    
else
    SpikeTimes = SpikeTimes(SpikeTimes >= trackingts(1) & SpikeTimes <= trackingts(end));
    session_duration =  trackingts(end) - trackingts(1);
end

%tossedspikes=length(session.units(u).ts)-length(SpikeTimes);
% if suppress<2
%     fprintf('tossing %d of a total %d spikes \n', tossedspikes, length(session.units(u).ts));
% end

averagerate=length(SpikeTimes)/session_duration;

end

