function [C,phi,S1,S2,S12] = calcSpikeFieldCoherence(spikeTimes,lfpData,trialTimes,params)
% function [C,phi,S1,S2,S12] = calcSpikeFieldCoherence(spikeTimes,lfpData,trialTimes,params)
%   This computes the spikeFieldCoherence by converting continuous data
%   into trialwise data.

%%%%%%%%%% *** YOU NEED CHRONUX FOR THIS TO WORK *** %%%%%%%%%%%%%

% prep these data:
%1. turn spike times into cell array locked to 0
%2. turn LFP into matrix with events
%3. input parameters
%4. get our data!


% start with LFP, its at a lower sampling frequency, so we want to snap
% triastarts to this.

% for each trial, find the first lfp ts that is after that start, pull that
% datapoint and the next n till trial end
% adjust that trialstart and end to those timestamps used
% move on.
lfpFS=mode(diff(lfpData(:,1))); % get most common sampling interval (which should just be the real one)
trialDur=ceil(namean(diff(trialTimes,1,2))/lfpFS); % need to make this matrix square!

lfpMat=nan(length(trialTimes),trialDur); % preallocate

for i=1:length(trialTimes(:,1))
    newStart=find(lfpData(:,1)>trialTimes(i),1,'first'); % find my start time
    lfpMat(:,i)=lfpData(newStart:newStart+trialDur,2); % fill in mat
    trialTimes(:,i)=[newStart newStart+(trialTimes(i,2)-trialTimes(i,1))]; % push trial window to the LFPts instead
end

% pull spike times for that trial matrix
[~,~,~,~,~,foundevspikes]=event_spikes(spikeTimes, trialTimes(:,1), 0, trialTimes(:,2));


%params=struct with: tapers, pad, Fs, fpass, err, trialave



[C,phi,S12,S1,S2,f,zerosp,confC,phistd,Cerr]=coherencycpt(lfpMat,foundevspikes,params)


end

