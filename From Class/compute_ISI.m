%Return the ISIs for the data d, not the autocorrellogram though**
%
%INPUT:
%  d = spike train data [n_trials, time].
%
%OUTPUT:
%  ISI = ISI values collapsed for all trials.
%
%MAK, Nov 2012.

function ISI = compute_ISI(d)

  n_trials = size(d,1);

  ISI = [];
  for k=1:n_trials
      spike_times = find(d(k,:) == 1);
      isi0 = diff(spike_times);
      ISI = [ISI, isi0];
  end
  
end