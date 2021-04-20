function [ripEvents,ripEEG] = GetRippleEvents(eeg,varargin)
% filters eeg data at ripple band frequency, then outputs ripple events
% INPUTS
%   REQUIRED
%   eeg: struct with fields:
%       data: raw data values (voltage)
%       samprate: should be 1000, or 1500
%       starttime: usually 0, or first increment (0.001 or 0.000667)
%       descript: 'recorded at... filtered at...'
%   OPTIONAL:
%   minZ: min z value for calculating duration and start-end of ripple, in
%       z, usually 1 or 1.5 (sd above mean)
%   minPeakZ:  minimum peak intensity for it to be a ripple, usually 3-4 (sd
%       above mean)
%   smoothSpan: Width of sd of smoothing kernel, generally 3-6 msec, for
%       smoothing the envelope of our ripple power.
%   minDelay: min off time before next ripple, if its below this, the
%       ripples are combined
% OUTPUTS
%   ripEvents: Struct iwth fields:
%       Area: the duration basically (i'll convert to msec
%       Centroid: The center index (look to LFPts)
%       PixelValues: Z values
%       MeanIntensity: mean Z score
%       BoundingBox: start and end indices (in LFPts)

p = inputParser;

addOptional(p,'minZ',1);
addOptional(p,'minPeakZ',3);
addOptional(p,'smoothSpan',5);
addOptional(p,'minDelay', 50);
parse(p);

minZ=p.Results.minZ;
minPeakZ=p.Results.minPeakZ;
smoothSpan=p.Results.smoothSpan;

try
    load('ripple1k.mat');
catch
    [myfile,mypath]=uigetfile('','Find the filter file');
    load(fullfile(mypath,myfile));
end

ripEEG=filtereeg2(eeg,ripple1k,'int16',0); % the guts here is filtfilt2(kernel,1,rawdata

ripraw=ripEEG.data(:,3); % 3 is amplitude envelope

ripsmooth=SmoothMat2(ripraw,[1, smoothspan*4],smoothspan); % smooth out 4x the std
ripZ=zscore(ripsmooth); % zscore our envelope
ripEEG.data(:,4)=ripZ;

ripEvents=regionprops(ripZ>putativepwr,ripZ,'Area',...
    'Centroid','PixelValues','MeanIntensity','BoundingBox');

% cut rips who dont peak high enough
ripEvents=ripEvents(cellfun(@(a) max(a)>peakpwr, {ripEvents.PixelValues}));

% now cut rips that are too short
ripEvents=ripEvents(cellfun(@(a) a>minRipLength*eeg.samprate, {ripEvents.Area}));

% now clean up the fields, 
% and if necessary cut or combine rips too soon after the previous ripple
%lastRipEnd=nan; killist=[];

for i=1:length(ripEvents)
    ripEvents(i).Centroid=ripEvents(i).Centroid(end);
    ripEvents(i).BoundingBox=[ripEvents(i).BoundingBox(2) ripEvents(i).BoundingBox(2)+ripEvents(i).BoundingBox(4)];
    ripEvents(i).Area=ripEvents(i).Area*eeg.samprate;
    
    % insert code here to see if the beginning of this rip is too close to
    % the end of the last rip
    %{
    if potrips(i).BoundingBox(1)<(lastRipEnd+minDelay)
        % combine
        potrips(i-1).BoundingBox(2)=potrips(i).BoundingBox(2);
    end
    lastRipEnd=
    %}
end

end

%% code base that inspired this function
%{

from supplement of Coordinated Excitation and Inhibition of Prefrontal Ensembles
 during Awake Hippocampal Sharp-Wave Ripple Events

All analyses were performed using custom code written in Matlab (MathWorks, Inc.).
SWRs were detected using LFPs filtered in the 150-250 Hz range on multiple CA1
tetrodes as previously described (Karlsson and Frank, 2009). Increases in power in the
ripple band were detected using a 3 s.d. criterion. A speed criterion of < 4 cm/sec was
also applied to avoid any spurious detection of high frequency fast gamma episodes
which occur during running (Kemere et al., 2013). Further, we discarded all SWRs that
occurred within 1 sec of a previously detected SWR, ensuring unambiguous assignment
of PFC spikes with respect to SWR onset time. A total of 19,888 SWRs were recorded
over 53 days (mean ± sem = 381 ± 18 per day). The average duration of SWRs was
96.7 ± 1.2 ms (n = 19,888 SWRs), and the average SWR rate was 0.12 ± 0.01 Hz (n =
53 days of recording).


from tang 2017, hippocampal-prefrontal reactivation during learning is
stronger in awake as compared to sleep states
SWR detection and modulation. SWRs were detected as previously de-
scribed (Cheng and Frank, 2008; Jadhav et al., 2016). Briefly, LFPs from
multiple CA1 tetrodes were filtered into the ripple band (150�250 Hz) 
and the envelope of bandpassed LFPs was determined by Hilbert transform. 
SWRs were initially detected when the envelope exceeded a threshold 
(mean ? 3 SD) on at least one tetrode. Detection ofSWRs was performed only 
when the animals� head speed was ?4 cm/s. SWR events were then defined as 
times around the initially detected events during which the envelope exceeded 
the mean. The duration of a SWR event was defined as the difference between 
its onset and offset and the amplitude was defined in terms of exceeded SDs 
above the mean (SWR properties summarized in Fig. 4F). For sleep analysis, 
only SWRs that occurred in SWS periods were included. For SWR-triggered 
rasters, peri- stimulus time histograms (PSTHs), and spectral analysis, 
only SWRs separated from others by at least 500 ms were included


% pull some random data in
lfpchan=1;
lfpchan=9;

eeg.data=SuperRat(45).LFP(lfpchan).data;

eeg.samprate=1000;
eeg.starttime=SuperRat(45).LFPts(1);
%eeg.data= []; % raw data in column vector
eeg.descript= 'data recorded at 1khz analog filtered below 300hz';
rippledata=filtereeg2(eeg,ripple1k,'int16',0); % the guts here is filtfilt2(kernel,1,rawdata

% then once you get the rippledata you need to convolve with a gaussian of
% 4 ms std, and then extract 4 std above mean patches using extractevents

ripraw=rippledata.data(:,3); % 3 is amplitude envelope
smoothspan=5; % msec
ripsmooth=SmoothMat2(ripraw,[1, smoothspan*4],smoothspan); % smooth out 4x the std

figure; plot(ripraw(1000:2000));
hold on;
plot(ripsmooth(1000:2000));


ripZ=zscore(ripsmooth);

% okay the other way to do this is using the image processing toolbox
% signal, thresh, mean, 0, min duration, another 0
%rips=extractevents(ripraw, 3, 0, 0, 50, 0);

% first grab the blobs where it passes threshold
putativepwr=1.5;
peakpwr=3;

% find all blobs above our size
potrips=regionprops(ripZ>putativepwr,ripZ,'Area',...
    'Centroid','PixelValues','MeanIntensity','BoundingBox');

% remove those without a large enough peak
potrips=potrips(cellfun(@(a) max(a)>peakpwr, {potrips.PixelValues}));

% now remove those that are too short
minRipLength=.04;
potrips=potrips(cellfun(@(a) a>minRipLength*eeg.samprate, {potrips.Area}));
% ripple event detection
% ripple events are first detected as events hwen the power exceeds 3-4 sd above the mean
% the duration is flood-fill until it goes to 1.5, 1, or 0 sd above ,mean
% so the way i will do it is grab ripples above 1sd above mean, and then
% trim those whose peak is below 3 sd above mean.

% now kill small ones
%potrips=potrips(potrips.test);
for i=1:length(potrips)
    potrips(i).Centroid=potrips(i).Centroid(end);
    potrips(i).BoundingBox=[potrips(i).BoundingBox(2) potrips(i).BoundingBox(2)+potrips(i).BoundingBox(4)];
    % insert code here to see if the beginning of this rip is too close to
    % the end of the last rip
end

%}