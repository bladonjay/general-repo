function beeps = spikeBeeps(Fs, n, duration)

freqs = linspace(1200,3200,n);
beeps = single(zeros(n,ceil(duration*Fs)));
ts = (0:ceil(duration*Fs)-1)/Fs;
for b = 1:n
    beeps(b,:) = cos(freqs(b)*pi*ts).*sin(ts*pi/duration);
end

end