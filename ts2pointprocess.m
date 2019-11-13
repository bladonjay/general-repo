function [spikevector,times] = ts2pointprocess(timestamps,varargin)
%function [spikevector,times] = ts2pointprocess(timestamps,varargin)
% this function takes a vector of timestamps and turns it into a vector of
% trues and falses and decimates your timestamps to a decent value
 
ts=timestamps;

% parse varargin
p=inputParser;

% factor to which spike times will be rounded/accurate (e.g. 10000 = 1/10000 of a second)
addOptional(p,'timebin',.001);

% Maximum lag to use when calculating autocorrelegram
addOptional(p,'max_lag_ms',200);

parse(p,varargin{:});


timebin=p.Results.timebin;
max_lag_ms = p.Results.max_lag_ms; 
%max is 10000 with only 2GB RAM available to MATLAB before crashing
% Nowadays we have like 8GB
max_num_spikes = 300000; 
 
max_lag_scaled = (1/timebin) * (max_lag_ms/1000);

% which timestamps to use
t_use = ceil(ts(1:min(max_num_spikes,length(ts)))*(1/timebin));
 
spikevector = zeros(1,max(t_use));

spikevector(t_use) = ones(1,length(t_use));

times=(timebin:timebin:max(t_use)*timebin);
end

