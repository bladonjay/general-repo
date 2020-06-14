%%  Compute the average envelope versus phase.  %%
%   From Buzsaki et al, Neuroscience, 2003.
%
%   INPUTS:
%    a = amplitude envelope signal.
%    p = phase signal (in radians).
%
%   OUTPUTS:
%    Plot the average and standard deviation of the amplitude envelope versus %    the phase
%
%   MAK.  Nov 12, 2007.
 
function average_envelope_versus_phase(a, p, varargin)
 
  p = p * 360.0 / (2.0*pi);                 %Use phase in degrees.
  
  a_mean = zeros(37,1);                     %Divide the phase into 37 bins,
  a_std = zeros(37,1);                      %Each of width 10 degrees.
  angle = zeros(37,1);
  j=1;
  
  for k=-180:10:180
      indices = find(p >= k & p < k+10);  %Find phase values in a 10 degree band.
      a_mean(j) = mean(a(indices));       %Compute mean amplitude at these phases.
      a_std(j) = std(a(indices));         %Compute std of amplitude at these phases.
      angle(j) = k+5;                     %Label the phase bin with the center phase.
      j=j+1;
  end
  
  %Plot the mean envelope versus phase.
  
  plot(angle, a_mean, 'k', 'LineWidth', 1);
  hold on
  plot(angle, a_mean+a_std, '--', 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
  plot(angle, a_mean-a_std, '--', 'Color', [0.8,0.8,0.8], 'LineWidth', 1)
  hold off
  set(gca, 'FontSize', 18);
  xlim([-180 180])
  xlabel('Theta phase');  ylabel('Gamma envelope');
 
end

