%%  Compute the two-dimensional modulation index and plot the result.  %%
%   From Canolty et al, Science, 2006.
%   INPUTS:
%    d      = The unfiltered data.
%    srate  = The sampling rate (e.g., 1000 Hz)
%
%   OUTPUS:
%    A plot of the two-dimensional modulation index in grayscale.
%    mod2d  = The two-dim modulation index.
%    flow   = The frequency axis for the lower phase frequencies.
%    fhigh  = The frequency axis for the higher amplitude frequencies.
%
%   DEPENDENCIES:
%    modulation_index.m
%    eegfilt.m (From EEGLAB Toolbox, http://www.sccn.ucsd.edu/eeglab/)
%
%   MAK.  Nov 12, 2007.
 
function [mod2d, flow, fhigh] = modulation_index_2d(d, srate)
 
  flow1 = (2:1:40);                                 %The phase frequencies.
  flow2 = flow1+1.0;                                %2-41 Hz with 1 Hz steps.
  
  fhigh1 = (5:5:200);                             %The amp frequencies.
  fhigh2 = fhigh1+5.0;                            %5-205 Hz with 5 Hz steps.
  
  mod2d = zeros(length(flow1), length(fhigh1));
  
  for i=1:length(flow1)
      theta=eegfilt(d,srate,flow1(i),flow2(i));  %Compute the low freq signal.
      theta=theta(srate:length(d)-srate-1);           %Drop the first and last second.
      phase = angle(hilbert(theta));                  %Compute the low freq phase.
      ['Loops remaining = ', num2str(length(flow1)-i)]
      for j=1:length(fhigh1)
        gamma=eegfilt(d,srate,fhigh1(j),fhigh2(j));%Compute the high freq signal.
        gamma=gamma(srate:length(d)-srate-1);         %Drop the first and last second.
        amp = abs(hilbert(gamma));                 %Compute the high freq amplitude.

        %Compute the modulation index.
        [m_norm_length] = modulation_index(amp, phase, srate);
        mod2d(i,j) = m_norm_length;
      end
  end
 
  %Plot the two-dimensional modulation index.
  
  flow = (flow1 + flow2) / 2;
  fhigh = (fhigh1 + fhigh2) / 2;
  imagesc(flow, fhigh, mod2d', [2 10]);  colorbar;
  axis xy
  set(gca, 'FontSize', 18);
  xlabel('Phase Frequency [Hz]');  ylabel('Envelope Frequency [Hz]');
  colormap(1-gray)
 
end