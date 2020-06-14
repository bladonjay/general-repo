%%  Plot the unfiltered and filtered data, amplitude envelope and phase.
%
%  INPUTS:
%   i_start  = start index
%   i_stop   = stop index
%   dt       = sampling interval
%   s        = unfiltered data.
%   theta    = theta filtered data.
%   gamma    = gamma filtered data.
%   phase    = phase of the theta signal.
%   amp      = amplitude envelope of gamma signal.
%
%  OUTPUTS:
%   Plot the unfiltered and filtered data, amplituded envelope, and phase.
%
%  MAK. Nov 12, 2007.
 
function plot_filtered_data(i_start, i_stop, dt, s, theta, gamma, phase, amp)
 
  t = (0:i_stop-i_start)*dt;
 
  plot(t, s(i_start:i_stop)+8.5, 'k', 'LineWidth',1)
  ylim([3 10])
  
  hold on
  plot(t, theta(i_start:i_stop)+6, 'Color', [0.8, 0.8, 0.8], 'LineWidth',1.0)
  plot(t, phase(i_start:i_stop)./pi+6, 'k', 'LineWidth',1)
 
  plot(t, gamma(i_start:i_stop)+4, 'Color', [0.8, 0.8, 0.8], 'LineWidth',1.0)
  plot(t, amp(i_start:i_stop)+4.1, 'k', 'LineWidth',1)
  hold off
  
  set(gca, 'Ytick', [], 'FontSize', 18);
  xlabel('Time [s]')
 
end
