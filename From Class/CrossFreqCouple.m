function [Phaseamp,Metric]=CrossFreqCouple(signal,dt,lowfreq,hifreq,plotter)
% function [Phaseamp,Metric]=CrossFreqCouple(signal,dt,lowfreq,hifreq,plotter)
%calculate coherence between low freq phase and high freq amplitude
% signal is vector of voltages, time is timepoints, lowfreq will be a
% window [low high] and high freq will be the same, [low high]
% first four inputs required, plot not required

if ~exist('plotter','var')
    plotter=0;
elseif isempty(plotter)
    plotter=0;
end


%   Plot the data, and check through visual inspection whether CFC occurs,

N = length(signal);  Fs = 1/dt;          %Define useful parameters.
t = (0:N-1)*dt;%Define a time axis.
if plotter
    figure
    subplot(2,1,1);
    plot(t,signal);
    xlabel('Time [s]'); ylabel('signal');
    subplot(2,1,2);
    periodogram(signal, [], [], 1/dt) % really useful to use
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Step 1:  Filter the data into low and high frequency bands.  To filter
%   the data, we'll use a second order Butterworth filter.  (Note that
%   there are many possible choices here!)

% 5-7 hz i think
deg=2;                              %Sets the filter order.
Wn = [lowfreq(1)*2/Fs,lowfreq(2)*2/Fs];               %Define the low frequency window of interest.
[B,A] = butter(deg,Wn,'bandpass');  %Apply the filter to isolate the band.
dlo = filtfilt(B,A,signal);

% 60-140 hz
Wn = [hifreq(1)*2/Fs,hifreq(2)*2/Fs];            %Define the high frequency window of interest.
[B,A] = butter(deg,Wn,'bandpass');  %Apply the filter to isolate the band.
dhi = filtfilt(B,A,signal);
if plotter
    figure
    plot(t,signal);
    hold on
    plot(t,dlo, 'r');
    plot(t,dhi, 'g');
    hold off
    xlabel('Time [s]');  title('Data (blue), Low-freq (red), High-freq (green)')
end
%   We might notice that funny effects happen at the edges of the data.  So,
%   let's ignore these edge effects due to the filtering and concentrate on
%   the center of the data.

dlo = dlo(N/4:end-N/4-1);  % this removes the ends because the function breaks
dhi = dhi(N/4:end-N/4-1);   % down at the ends of the signal
taxis = t(N/4:end-N/4-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Step 2:  Compute the phase of the low frequency signal and the
%   amplitude envelope of the high frequency signal.  To do so, we use the
%   Hilbert transform.

phi = angle(hilbert(dlo));    %Phase of low frequency signal.
amp = abs(hilbert(dhi));      %Amplitude envelope of high frequency signal.


%   Let's check that these results look okay,
if plotter
    figure;  clf()
    subplot(2,1,1)
    plot(taxis, dlo)
    hold on
    plot(taxis, phi, 'g')
    hold off
    axis tight
    xlabel('Time [s]');  title('Low-freq and phase');
    
    subplot(2,1,2)
    plot(taxis, dhi)
    hold on
    plot(taxis, amp, 'r')
    hold off
    axis tight
    xlabel('Time [s]');  title('High-freq and amplitude envelope')
end
%%%%%%%%%%%%%%%%
%   Step 4 now we compare the two, and find the phase where the max
%   amplitude is, and we find difference between min and max
%%%%%%%%%%%%%
p_bins = (-pi:0.1:pi);              %Define the phase bins.

a_mean = zeros(length(p_bins)-1,1); %Vector to hold avg amp env results.
p_mean = zeros(length(p_bins)-1,1); %Vector to hold center of phase bins.

% this is the part that creates the histogram
for k=1:length(p_bins)-1
    pL = p_bins(k);                         %Phase lower limit for this bin.
    pR = p_bins(k+1);                       %Phase upper limit for this bin.
    indices = find(phi >= pL & phi < pR);   %Find phase values in this range.
    a_mean(k) = mean(amp(indices));         %Compute mean amplitude at these phases.
    p_mean(k) = mean([pL, pR]);             %Label the phase bin with the center phase.
end
% heres the metric we use, and the one we will BOOTSTRAP later
Metric = max(a_mean)-min(a_mean);                %The difference between the max and min modulation.

if plotter
    figure
    plot(p_mean, smooth(a_mean), 'ko-', 'LineWidth', 1);
    axis tight
    xlabel('Low frequency phase');  ylabel('High frequency envelope height difference');
    title(['Metric h=' num2str(Metric)])
end
    Phaseamp=[p_mean a_mean];

end
