function [filtered,filterparams] = LFP_Filter(samples,timestamps,varargin)

%{
Filter - Filter samples.

 USAGE

   filtered = bz_Filter(samples,<options>)

INPUT
   samples        samples given as a [t x channels] timeseries vector
   timestamps     either a t timestamps timeseries vector or a single value of sampling interval             
   <options>      optional list of property-value pairs (see table below)

   =========================================================================
    Properties    Values
   -------------------------------------------------------------------------
    'passband'    pass frequency range. [0 X] for low-pass, [X inf] for highpass
    'stopband'    stop frequency range
    'order'       filter order (number of cycles, default = 4)
    'attenuation'      filter attenuation (default = 20 db of filtering)
    'filttype'      choose filter type between 'cheby2' (default) and 'fir1'
    'nyquist'     nyquist frequency (default = 625), calculated 
                  automatically from samples.samplingRate for BUZCODE
    'FMAlegacy'   true/false, uses FMA legacy input style
                  (samples given as a list of (t,v1,v2,v3...) tuples)
    'fast'        true/false, uses FiltFiltM if true (doesn't work w/ new
                  version of matlab)
    'channels'    if input is a buzcode lfp structure with field
                  samples.channels, will only filter the selected
                  channels
    'intervals'   only filter in given intervals (only works for buzcode input)
   =========================================================================

OUTPUT
  filtered        -if the input is a timeseries vector, output is as well
                  -if the input is a buzcode structure, output is a
                   buzcode structure with the following fields:
      .data       the filtered data
      .phase      phase, calculated by the hilbert method
      .amp        amplitude, calculated by the hilbert method
      .timestamps     (s)
       .samplingRate   (Hz)
       .filterparms    a structure that holds the parameters used for
                       filtering, for future reference

 Copyright (C) 2004-2011 by Michaël Zugaro
 updated 2017 DLevenstein for buzcode
updated again 2020 john bladon personal code
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%}
%% Params and Defaults and outputs 
filtered=[];
filterprops=struct('passband',[]);

% parse inputs
p=inputParser;
p.addParameter('passband',[6 12], @(a) erfpassband);
p.addParameter('stopband',[],@(a) erfstopband);
p.addParameter('filttype','fir1',@(a) erffilttype);
p.addParameter('order',5,@(a) isnumeric(a));
p.addParameter('attenuation',20,@(a) isnumeric(a));
% Parse parameter list

parse(p,varargin{:});
			
filterprops=p.Results;
order=filterprops.order;
passband=filterprops.passband;
stopband=filterprops.stopband;
attenuation=filterprops.attenuation;
if length(timestamps) > 1
    samplingrate=diff(timestamps);
    nyquist = 0.5/median(samplingrate);
else
    nyquist = 0.5/timestamps;
end

% build filter
switch(filterprops.filttype)
    case 'cheby2'
        if ~isempty(passband)
            if passband(1) == 0
                [b a] = cheby2(order,attenuation,passband(2)/nyquist,'low');
            elseif passband(2) == inf
                [b a] = cheby2(order,attenuation,passband(1)/nyquist,'high');
            else
                [b a] = cheby2(order,attenuation,passband/nyquist);
            end
        else
            [b a] = cheby2(order,attenuation,stopband/nyquist,'stop');
        end
        warning('Cheby2 is often numerically unstable - if you get NaNs, use the ''fir1'' argument to bz_Filter');
    case 'fir1'
        %Order input to fir1 needs to be in samples, not cycles
        if ~isempty(passband)
            if passband(1) == 0
                filt_order = round(order*2*nyquist./passband(2));
                [b a] = fir1(filt_order,passband(2)/nyquist,'low');
            elseif passband(2) == inf
                filt_order = round(order*2*nyquist./passband(1));
                [b a] = fir1(filt_order,passband(1)/nyquist,'high');
            else
                filt_order = round(order*2*nyquist./passband(1));
                [b a] = fir1(filt_order,passband/nyquist);
            end
        else
            filt_order = round(order*2*nyquist./stopband(1));
            [b a] = fir1(filt_order,stopband/nyquist,'stop');
        end
    case 'butter'
        if ~isempty(passband)
            [b a] = butter(order,[passband(1)/nyquist passband(2)/nyquist],'bandpass');
        else
            [b a] = butter(order,stopband(1)/nyquist,'stop');
        end
end


% now filter the data
 for i = 1:size(data,2) % now for each col get the filtered data
        if ~fast
           filtered(:,i) = filtfilt(b,a,double(samples(:,i)));
        else
           filtered(:,i) = FiltFiltM(b,a,double(samples(:,i))); 
        end
 end


end

%{
add and remove overhang?
%Restrict to intervals, with overhang to remove edge effects at transitions
    %(then remove later)
    overhang = (order)./passband(1);
    overint = bsxfun(@(X,Y) X+Y,intervals,overhang.*[-1 1]);
    keepIDX = InIntervals(samples.timestamps,overint);
    samples.data = samples.data(keepIDX,:);
    filtered.timestamps = samples.timestamps(keepIDX);


 %Remove the overhang from intervals
    keepIDX = InIntervals(filtered.timestamps,intervals);
    filtered.data = filtered.data(keepIDX,:);
    filtered.hilb = filtered.hilb(keepIDX,:);
    filtered.amp = filtered.amp(keepIDX,:);
    filtered.phase = filtered.phase(keepIDX,:);
    filtered.timestamps = filtered.timestamps(keepIDX);

%}
%% all the error functions

function [rez]=erfpassband(band)
rez=(isnumeric(band)) && (length(band)==2) && (diff(band)>0);
if ~rez
    error('Incorrect value for ''passband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end
end

function [rez]=erfstopband(band)
rez=(isnumeric(band)) && (length(band)==2) && (diff(band)>0);
if ~rez
    error('Incorrect value for ''stopband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end
end

function [rez]=erffilttype(myfilt)
templates={'cheby2','fir1','butter'};
okfilt=ismember(lower(myfilt),lower(templates));
if any(okfilt)
    fprintf('using %s filter \n',templates{okfilt(find(okfilt,1,'first'))});
end
rez=any(okfilt);
end
    