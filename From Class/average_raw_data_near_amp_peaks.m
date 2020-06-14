%%  Compute the event-related averaged triggered by peaks in the gamma signal.
%   From Bragin et al, J Neurosci, 1995.
%
%   INPUTS:
%    gamma  = the high-frequency data (e.g., data filtered at gamma).
%    s      = the unfiltered signal.
%
%   OUTPUTS:
%    Plot the unfiltered data average around the gamma peaks (+/- 100 ms).
%    avg    = the unfiltered data located +/- 100 ms around each event detected.
%    lags   = the time axis around each event (+/- 100 ms).
%
%   DEPENDENCIES:
%    findpeaks.m (New in MATLAB 7.5 R2007b)
%
%   MAK.  Nov 14, 2007.
 
function [avg, lags] = average_raw_data_near_amp_peaks(gamma, s)
 
  win = 100;
 
  [pks, loc] = findpeaks(gamma, 'minpeakdistance', 100);
  
  iloc = find(loc > win & loc < length(s)-win);
  lags = (-win:win);
  loc = loc(iloc);
  
  avg = zeros(length(loc), 2*win+1);
  for i=1:length(loc)
      avg(i,:) = s(loc(i)-win:loc(i)+win);
  end
  
  plot(lags, sum(avg,1)/length(loc), 'k');
  axis tight
  ylim([-1.5, 1.5])
  set(gca, 'FontSize', 18);
  xlabel('Time [ms]');  ylabel('Average []');
  
end
