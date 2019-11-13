function spkmat = spikesMat(data, neurons, ts)
% Create an m x n matrix of number of spikes in each time bin.
% m is the length of ts. n is the length of neurons.

nts = numel(ts);
nneuron = numel(neurons);

spks = cell(nneuron,1);
if isfield(data.neurons(1),'timestamps')
    for ii = 1:nneuron
        spks{ii}(:,1) = data.neurons{neurons(ii)}.timestamps;
        spks{ii}(:,2) = ii;
    end
elseif isfield(data.neurons(1),'ts')
    for ii = 1:nneuron
        spks{ii}(:,1) = data.neurons(neurons(ii)).ts;
        spks{ii}(:,2) = ii;
    end
end

spks2 = cell2mat(spks);

% construct bins
dt = diff(ts);
bins = [ts(1)-dt(1)/2; ts(1:end-1)+dt/2; ts(end)+dt(end)/2];
bins(end) = bins(end)+eps(bins(end));

% count spikes in each bin
[~,spks2(:,1)] = histc(spks2(:,1),bins);
spks2(spks2(:,1)==0,:)=[];
spkmat = accumarray(spks2,1,[nts nneuron]);

end