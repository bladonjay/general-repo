function [faxis,Sxx,fNQ,T] = PlotPowerSpectrum(data,time,varargin)
%PlotPowerSpectrum plots the power spectrum of the data.  It uses a Hanning
%filter unless otherwise noted, and plots both decibels in a log frequency
%axis and a power spectrum in a linear x axis.  It also allows you to
%designate the freq window for the plot.

% INPUTS
% data should be a membrane potential
% time should start at 0 and be monotonically increasing

% OUTPUTS
% Faxis= x axis for power spectrum
% Sxx is the power
% fNQ is the nyquist freq( highest freq you can see)
% T= total time


p=inputParser;
addOptional(p,'hanning',1);
addOptional(p,'xlimit',[1 100]);
addOptional(p,'suppress',0);
parse(p,varargin{:});

xlimit=p.Results.xlimit;
hanning=p.Results.hanning;
suppress=p.Results.suppress;
t=time;
d=data(:);


% calculate some base parameters
T=max(t)-min(t);
df=1/T;   %Freq resolution

dt=t(2)-t(1);
fNQ=1/dt/2; % nyquist frequency


if hanning==1;
    d=d.*hann(length(d));
    titlename='hanning filtered power spectrum';
else
    titlename='square filtered power spectrum';
end

if ~(suppress)
    figure;
    plot(t,d);
end


xf = fft(d);                     %Define the frequency resolution.                      %Define the Nyquist frequency.
faxis = (df:df:fNQ);                 %Now, define the frequency axis,
Sxx = 2*dt^2/T * (xf.*conj(xf));    %Compute the power spectrum.
Sxx = Sxx(1:(length(Sxx)/2)-1);       %... and keep only power at pos freq.

axisend=min([length(faxis), length(Sxx)]);
faxis=faxis(1:axisend);
 Sxx=Sxx(1:axisend,:);
if ~(suppress)
    figure;
    subplot(1,2,1);
    semilogx(faxis, Sxx)     % 10*log(Sxx)   %Plot the power versus frequency,
    xlabel('Frequency (Hz)');           %... and label the axes.
    ylabel('Power in db');
    title(titlename)
    xlim(xlimit)
    
    subplot(1,2,2);
    plot(faxis, Sxx)        %Plot the power versus frequency,
    xlabel('Frequency (Hz)');           %... and label the axes.
    ylabel('Power');
    title(titlename);
    xlim(xlimit);
end

end

