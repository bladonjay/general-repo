function [cor, lags] = CrossCorr(ts1, varargin)
% [acor, lag] = CMBHOME.Spike.CrossCorr(ts1)
% [xcor, lag] = CMBHOME.Spike.CrossCorr(ts1, ts2, varargin)
%
% Calculates cross correlation (or auto correlation) for vectors ts1 and
% ts2 of spike times. ts1 and ts2 are not binned, or binary, spike trains but
% rather the time that spikes occured.
%
% If ts1 is the only argument, the autocorrelation is produced (also if
% ts1==ts2). In this case, all zeros in the set of latencies will be removed
% before calculating the xcorr histogram.
%
% ARGUMENTS
%   ts1         vector of spike times (seconds)
%   ts2         (optional) vector of spike times
%
% PROPERTIES
%   'binsize'       (default=.01 seconds) binsize for cross correlation
%   'norm'          (default='prob') determines normalization of cross
%                   correlogram. Can be 'prob', 'count', 'unbiased'. To use
%                   the unbiased option, 'epoch' property must be passed.
%   'lag'           (default [-.4 4] seconds) included to speed up algorithm.
%                   this defines the upper and lower limit in the peristumulus
%                   time histogram. can be multiple rows.
%   'suppress_plot' (default 1) 1 with plot, 0 without plot
%   'epoch'         2 element vector indicating start and stop times of recording window
%
% RETURNS
%   cor         a col. vector of cross correlation values corresponding to 'lag'
%   lag         a col. vector of center-aligned lag values (seconds)
%
% alex wiltschko and andrew bogaard

p = inputParser;

p.addRequired('ts1');
p.addOptional('ts2', NaN, @(c) isnumeric(c));
p.addParamValue('norm', 'prob', @(c) ischar(c));
p.addParamValue('binsize', .01, @(c) all(c>0));
p.addParamValue('lag', [-.4 .4], @(c) size(c,2)==2 && all(c(:,2)-c(:,1))>0);
p.addParamValue('suppress_plot', 1, @(c) numel(c)==1);
p.addParamValue('epoch', [], @(c) numel(c)==2);

p.parse(ts1, varargin{:});

ts1 = p.Results.ts1(:); % make col vectors
ts2 = p.Results.ts2(:);
norm = p.Results.norm;
binsize = p.Results.binsize;
lag = p.Results.lag;
suppress_plot = p.Results.suppress_plot;
epoch = p.Results.epoch;

ac = 0;

if isnan(ts2), ts2 = ts1; end

if isequal(ts1,ts2), ac = 1; end % is autocorr

if length(binsize) ~= size(lag,1), error('lag and binsize must be NxM and N sizes, respectively'); end

if isempty(ts2) || isempty(ts1),
    warning('one of the spikes trains is empty. returning zero xcorr.');
    lags = lag(1)+binsize/2:binsize:lag(2)-binsize/2;
    lags = lags(:);
    cor = zeros(size(lags));
    return
end

db = nan(length(ts1), 3);

psth = [];

s1 = 1;

spkind = 1;

l1 = min(lag(:));
l2 = max(lag(:));

while spkind <= length(ts1)
   
    s = s1;
    
    while ts2(s) < ts1(spkind)+l1 && s < length(ts2)
    
        s = s+1;
    
    end
    
    s1 = s;
    
    f = s;
    
    while ts2(f) <= ts1(spkind)+l2
        
        f = f+1;
        
        if f>length(ts2), break; end
        
    end
        
    if ts2(s)<=ts1(spkind)+l2
        db(spkind, :) = [s f-1 ts1(spkind)];
    end
    
    spkind = spkind+1;
    
end

dspk = diff(db(:,1:2), 1, 2);

for i = 0:max(dspk)
    
    where = dspk>=i;
    
    tf = db(where,1)+i;
    
    psth = [psth; ts2(tf)-db(where,3)];

end

if ac, psth(psth==0) = []; end % remove zeros in autocorrelation

cor = cell(size(lag, 1), 1);
lags = cell(size(lag, 1), 1);

for i = 1:size(lag,1)

    [cor{i}, lags{i}] = HistAndPlot(lag(i,:), binsize(i), psth, norm, ts1, ts2, suppress_plot);
    
end

if length(cor)==1 % return numerical array if only one element long
    
    cor = cor{i};
    lags = lags{i};
    
end

function [cor, lags] = HistAndPlot(lag, binsize, psth, norm, ts1, ts2, suppress_plot)

lags = lag(1)+binsize/2:binsize:lag(2)-binsize/2;
cor = hist(psth, lags);

if strcmp(norm, 'unbiased')
   
    if isempty(epoch)
        
        disp('''epoch'' must be defined for unbiased normalization, no normalization performed')
    
    else
        
        L = floor(diff(epoch)/binsize);

        normc = abs(lags)/binsize;
        
        cor = cor./(L-normc);
        
    end
    
elseif strcmp(norm, 'prob')
    
    cor = cor / min(length(ts1), length(ts2));
    
end % otherwise, it will be count

if ~suppress_plot
    bar(lags*1000, cor, 1, 'k'), xlim(1000*[lag(1) lag(2)]), hold on
    line([0 0], ylim, 'linestyle', ':', 'color', [.7 .7 .7]);
    set(gca, 'fontsize', 8, 'box', 'off');
end

cor = cor(:);

lags = lags(:);

