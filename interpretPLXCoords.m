function [c, vtmodeout, vtmode] = interpretPLXCoords(ts, sv, vtmodein, forcemode)
% interpretPLXCoords - Interpret CinePlex video tracking data
%
% [c, vtmode] = interpretPLXCoords([ts sv])
% [c, vtmode] = interpretPLXCoords(ts, sv)
%
% Input:
%   ts - array of timestamps (in clock ticks)
%   sv - array of strobed event values (produced by CinePlex)
%   tvmodein - the anticipated vtmode of the input data
%   forcemode - force the output mode to match the anticipated mode
%
% Output:
%   c      - The coordinates
%            c(:, 1) - timestamp
%            c(:, 2) - x1
%            c(:, 3) - y1
%            c(:, 4) - x2 or motion (if present)
%            c(:, 5) - y2 (if present)
%            c(:, 6) - x3 (if present)
%            c(:, 7) - y3 (if present)%
%
%   vtmodeout - The video tracking mode used for the output coordinates.
%   vtmode - The video tracking mode determined based on the input data.
%            0 = UNKNOWN
%            1 = CENTROID        - 1 set of coordinates, no motion
%            2 = CENTROID/MOTION - 1 set of coordinates, with motion
%            3 = LED_1           - 1 set of coordinates
%            4 = LED_2           - 1 set of coordinates
%            5 = LED_3           - 1 set of coordinates
%            6 = LED_12          - 2 sets of coordinates
%            7 = LED_13          - 2 sets of coordinates
%            8 = LED_23          - 2 sets of coordinates
%            9 = LED_123         - 3 sets of coordinates
%
% AUTHOR: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% Copyright (c) 2011, Benjamin Kraus
% $Id: interpretPLXCoords.m 3412 2012-01-18 23:37:07Z bkraus $

% Check if both timestamps and strobbed values were provided in the first
% input argument.
if((nargin == 1 || isempty(sv)) && size(ts,2) == 2);
    sv = ts(:,2); ts = ts(:,1);
elseif(nargin >= 2 && ~all(size(sv)==size(ts)))
    error('InterpretPLXCoords:ArraySizeMismatch',...
        '''ts'' and ''sv'' inputs must be same size.');
elseif(nargin < 2 && ~isempty(ts))
    error('InterpretPLXCoords:MissingValues',...
        'Cannot find strobbed values.');
elseif(isempty(ts)); sv = [];
end

% Convert strobed values into uint16 for easier manipulation later.
sv = uint16(sv);

% Determine the number of values transfered.
n = size(ts,1);

% WARNING: ****MATLAB uses 1 based numbering for bits****
% The notes below will differ from the CinePlex documentation because the
% CinePlex documentation uses 0 based numbering for bits.

% Bit 15 is always 1, it not, something is wrong.
check = bitget(sv,15);

if(any(check~=1));
    error('InterpretPLXCoords:InvalidData','Check bit (14) is not set.');
end

% Bit 14 to 11 determine the source of the data value.
try
    code = bitshift(sv,-10,4);
catch
    warning('bitshift didnt work so trying it on 16 bit output');
    code = bitshift(sv,-10,'uint16');
end
% 0000 ( 0) - ContrastX
% 0001 ( 1) - ContrastY
% 0100 ( 4) - Motion
% 0101 ( 5) - Color1X
% 0110 ( 6) - Color1Y
% 0111 ( 7) - Color2X
% 1000 ( 8) - Color2Y
% 1001 ( 9) - Color3X
% 1010 (10) - Color3Y

% Bit 9 to 1 store the actual value.
val = bitand(uint16(sv),1023);

% Determine whether the timestamps are in seconds or in ticks.
% If there is only 1 (or no) timestamps, then these numbers don't matter.
if(n == 0); df = []; rate = 0;
elseif(n == 1); df = inf; rate = 0;
    
% If there are more than two values, then use 'mode' to determine the most
% common time step, and assume this is the transmission rate of the
% coordinates. Any time step larger than this is considered a new timestamp
else df = [inf; diff(ts)]; rate = mode(df);
end

% The cutoff is the maximum time gap between sets of coordinates.
cutoff = rate*2;

% Initialize the output matrix
c = nan(sum(df>cutoff), 12);

% Copy values into the output matrix based on coordinate type.
cur = 0;
for ii = 1:n
    if(df(ii)>cutoff); cur = cur + 1; c(cur,1) = ts(ii); end
    c(cur,code(ii)+2) = val(ii);
end

% Determine which sources are present and determine the video tracking mode
sources = any(~isnan(c),1);
sources = sources(1,2:end);

if(sources(3) || sources(4))
    vtmode = 0; % Undefined codes, should be empty.
elseif(sources(1) && sources(2)) % ContrastX and ContrastY
    if(~any(sources(5:11)))
        vtmode = 1; % Mode 1 = Contrast X and Y without Motion
    elseif(sources(5) && ~any(sources(6:11)))
        vtmode = 2; % Mode 2 = Contrast X and Y with Motion
    else vtmode = 0; % Mode 0 = Unknown
    end
elseif(sources(1) || sources(2) || sources(5))
    vtmode = 0; % Contrast X or Y without both (invalid)
elseif(all(sources([6 7          ])) && ~any(sources([    8 9 10 11])))
    vtmode = 3; % Mode 3 = LED 1
elseif(all(sources([    8 9      ])) && ~any(sources([6 7     10 11])))
    vtmode = 4; % Mode 4 = LED 2
elseif(all(sources([        10 11])) && ~any(sources([6 7 8 9      ])))
    vtmode = 5; % Mode 5 = LED 3
elseif(all(sources([6 7 8 9      ])) && ~any(sources([        10 11])))
    vtmode = 6; % Mode 6 = LED 1 and 2
elseif(all(sources([6 7     10 11])) && ~any(sources([    8 9      ])))
    vtmode = 7; % Mode 7 = LED 1 and 3
elseif(all(sources([    8 9 10 11])) && ~any(sources([6 7          ])))
    vtmode = 8; % Mode 8 = LED 2 and 3
elseif(all(sources(6:11)))
    vtmode = 9; % Mode 9 = LED 1, 2, and 3
else vtmode = 0; % Mode 0 = Unknown/Invalid
end

% Translate 'vtmodein' to a logical 'sources' array.
if(exist('vtmodein','var')==1 && ~isempty(vtmodein))
    switch vtmodein
        case 0; sourcesin = false(1,11);                      % Mode 0 = Unknown
        case 1; sourcesin = logical([1 1 0 0 0 0 0 0 0 0 0]); % Mode 1 = Contrast X and Y without Motion
        case 2; sourcesin = logical([1 1 0 0 1 0 0 0 0 0 0]); % Mode 2 = Contrast X and Y with Motion
        case 3; sourcesin = logical([0 0 0 0 0 1 1 0 0 0 0]); % Mode 3 = LED 1
        case 4; sourcesin = logical([0 0 0 0 0 0 0 1 1 0 0]); % Mode 4 = LED 2
        case 5; sourcesin = logical([0 0 0 0 0 0 0 0 0 1 1]); % Mode 5 = LED 3
        case 6; sourcesin = logical([0 0 0 0 0 1 1 1 1 0 0]); % Mode 6 = LED 1 and 2
        case 7; sourcesin = logical([0 0 0 0 0 1 1 0 0 1 1]); % Mode 7 = LED 1 and 3
        case 8; sourcesin = logical([0 0 0 0 0 0 0 1 1 1 1]); % Mode 8 = LED 2 and 3
        case 9; sourcesin = logical([0 0 0 0 0 1 1 1 1 1 1]); % Mode 9 = LED 1, 2, and 3
        otherwise; sourcesin = false(1,11);
            warning('InterpretPLXCoords:InvalidVTmode','Ignoring invalid input VTmode');
    end
else
    vtmodein = vtmode;
    sourcesin = sources;
end

% Compare the actual 'vtmode' with the anticipated 'vtmodein'
% If coordinates were anticipated but not found, and no coordinates were
% found that were not anticipated, then use 'vtmodein' instead of 'vtmode'.
if(exist('forcemode','var')==1 && ~isempty(forcemode) && forcemode)
    vtmodeout = vtmodein;
    sources = sourcesin;
elseif(any(sourcesin & ~sources) && ~any(sources & ~sourcesin))
    vtmodeout = vtmodein;
    sources = sourcesin;
else
    vtmodeout = vtmode;
end

% Keep only those entries of 'c' that are specified by the tracking mode.
c = c(:,[true sources]);

end
