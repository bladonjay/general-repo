
%%  Correlate the amplitude envelope with the signal and plot results.
%   From Hentschke et al, Eur J Neurosci, 2007
%
%   INPUTS:
%    signal         = the theta filtered signal.
%    gamma_envelope = the envelope of the gamma filtered signal
%
%   OUTPUTS:
%    Plot the CC of the raw data versus the CC of the surrogate data.
%    z-score    = the z-score of the max CC of the raw data versus the surrogates.
%    pval       = the corresponding p-value.
%
%   DEPENDENCIES:
%    shuffle_envelope_signal.m
%
%   MAK.  Nov 12, 2007.
 
function [zscore, pval] = correlate_envelope_signal(signal, gamma_envelope)
 
  [C_raw, lags] = xcorr(signal, gamma_envelope, 500, 'coeff');%Compute CC of raw data.
  
  nsurrogate = 100;                                           %Use 100 surrogates.
  C_surrogate = zeros(nsurrogate, length(C_raw));
  mx = zeros(nsurrogate,1);
  
  for k=1:nsurrogate                                          %For each surrogate,
      sf = shuffle_envelope_signal(gamma_envelope);           %Shuffle the envelope.
      C_surrogate(k,:) = xcorr(signal, sf, 500, 'coeff');     %Compute the CC of surr.
      mx(k) = max(abs(C_surrogate(k,:)));                       %Find the peak CC.
  end
  
  mx_raw = max(abs(C_raw));                                 %Find peak CC of raw data.
  [m_surrogate, std_surrogate] = normfit(mx);               %Fit surr. data.
  [h pval ci zscore] = ztest(mx_raw, m_surrogate, std_surrogate);  %Compute z-score.
 
  %Plot the raw CC versus the ensemble of surrogate of CCs.
  
  plot(lags, C_raw, 'k')
  hold on
  for i=1:100
      plot(lags, C_surrogate(i,:), 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.5')
  end
  plot(lags, C_raw, 'k')
  hold off
  set(gca, 'FontSize', 18);
  ylim([-0.55 0.55])
  xlabel('Time [ms]')
  ylabel('Cross Correlation')
  
end
