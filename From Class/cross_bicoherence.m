%Compute the bicoherence of data matrix d.
%The format of d is . . . d = d[number of trials, time of trial]
%f0 = The sampling frequency.
%fmax = the maximum frequency of interest.
%
%June 20, 2007.  MAK


% uses algorithm

%   (|<Xj1, Xj2, conj(Yj1+j2)>|^2) / (<Xj1>^2 <Xj2>^2 <Yj1+j2>^2))



function [b, faxis] = cross_bicoherence(d, d2, f0, fmax)

  sz = size(d);  %Determine the size to define useful parameters.
  K = sz(1);     %Number of trials.
  N = sz(2);     %Number of indices per trial.

  faxis = (0:N/2) / (N)*f0;       %Frequency axis, probably in Hz.
  good = find(faxis < fmax);      %The part of the frequency axis to consider.
  
  b = zeros(length(good), length(good));   %Variable to hold the results.
  faxis = faxis(good);                     %Keep only faxis of interest.
  
  % hold space for numerator
  numerator = complex(zeros(length(good)), zeros(length(good)));
  pow = zeros(N,1);
  pow1= zeros(N,1)
  
  % for k number of trials
  for k=1:K
      
      x = fft(hann(N).*d(k,:)'); %Take the FFT of segment with Hanning.
      y = fft(hann(N).*d2(k,:)');
      numTemp = complex(zeros(length(good)), zeros(length(good)));

      for f1=1:length(good)         %For each freq pair, compute numerator.
          for f2=1:length(good)
              numTemp(f1,f2) = x(f1)*x(f2)*conj(y(f1+f2));
          end
      end
      
      numerator = numerator + numTemp / K;   %Compute trial average of numerator.
      pow = pow + x.*conj(x) / K; %Compute FFT squared, and avg over trials.
      pow1= pow1 + y.*conj(y) / K; % do same for y
  end

  for f1=1:length(good)                      %Compute bicoherence.
      for f2=1:length(good)
          b(f1,f2) = abs(numerator(f1,f2)) / sqrt(pow(f1)*pow(f2)*pow1(f1+f2));
      end
  end

end  