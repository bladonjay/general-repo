%%  Compute the average power spectrum and plot the results.  %%
%
%   INPUTS:
%  	s    = signal of interest.
%  	dt   = sampling interval.
%
%   OUTPUTS:
%  	Plot the power spectrum averaged over 1000 ms non-overlapping windows.
%
%  MAK.  Nov 12, 2007.
 
function compute_and_plot_psd(s, dt)
 
  window = 1000;
  nensemble = floor(length(s)/window);
 
  avgPSD=0;
  for k=1:nensemble
    dat = s((k-1)*window+1:k*window);
    power = abs(fft(dat.*hann(length(dat))'-mean(dat.*hann(length(dat))')));
    power = power(1:length(dat)/2+1);
    power = 10.0*log10(power / max(power));
    avgPSD = avgPSD + power;
  end
  df = 1.0 / (length(dat)*dt);
  faxis = (0:length(dat)/2)*df;
  plot(faxis, power, 'k', 'LineWidth', 1);  xlim([df, 200.0]);  ylim([-25 0])
  set(gca, 'FontSize', 18);
  xlabel('Freq [Hz]');  ylabel('Power [dB]');
 
end
