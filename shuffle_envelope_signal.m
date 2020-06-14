%shuffle_envelope_signal.m

%% Implement the shuffling procedue from Hentschke et al, 2007.
%
%   INPUTS:
%    s = signal to shuffle.
%
%   OUTPUTS:
%    sf = the shuffled signal. %
%  MAK.  Nov 12, 2007.
 
function sf = shuffle_envelope_signal(s)
 
  ran = randn(length(s),1);             %Create Gaussian signal.
  [dumb, iran] = sort(ran);             %Rank order the gaussian.
  s0 = s(iran);                         %Sort s by rank of gaussian.
  ffts0 = fft(s0);                      %Compute the FFT of s0.
  ffts = fft(s);                        %And the FFT of s.
  As = abs(ffts);                       %Use the amplitude of s.
  phis0 = angle(ffts0);                 %And the angle of s0.
  ffts0new = As.*exp(i*phis0);          %To create new signal in freq domain.
  s1 = ifft(ffts0new);                  %Take inverse FFT to return to time.
  [dumb, is] = sort(s);                 %Rank order the ORIG signal.
  [dumb, is1] = sort(s1);               %Rank order the NEW signal.
  s1(is1) = s(is);                      %Replace rank ordered amp of new
                                        %with rank ordered amp of orig.
  sf = s1;
 
end


