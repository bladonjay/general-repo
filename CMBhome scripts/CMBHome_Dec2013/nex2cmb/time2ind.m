function inds = time2ind(x, xi)
%TIME2IND Finds the indices that correspond to a set of timestamps.
%   TIME2IND returns an M x P array of indices that correspond to the
%   M x P array of time values given in the value XI, based on the
%   reference list of timestamps provided in the value X. Unlike
%   timerange2ind, the index of the closest timestamp in X is used,
%   regardless of whether that timestamps is just before or just after the
%   timestamp in XI.
%
%   INDS = TIME2IND(X, XI)
%   
%   Input:
%       X: An M x N array where the first column represents an M x 1 list
%       of timestamps.
%       XI: An M x P array of timestamps.
%
%   Output:
%       INDS: An M x P array of indices. Each index represents the location
%       of the closest timestamp in the array X of the corresponding
%       timestamp in XI.
%
%   See also

% Written by Benjamin Kraus (bkraus@bu.edu) 2007

    lengthx = length(x(:,1));
    inds = interp1(x(:,1), (1:lengthx)', xi, 'nearest');

    inds(xi < x(1,1)) = 0;
    inds(xi > x(end,1)) = 0;
end