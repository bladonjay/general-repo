function [coords, mode] = readPLXCoords(filename)

if(nargin == 0)
    [filename,pathname] = uigetfile({'*.plx','PLX-files'},...
        'Select PLX file to load...');
    filename = fullfile(pathname,filename);
end

assert(ischar(filename) && exist(filename,'file')==2,...
    [mfilename ':InvalidFilename'],'Invalid filename or non-existent file.');

% Find the position data stored in the PLX file. It is assuming that the
% data came from a Plexon system, which stores the position data as the
% strobbed event number 257.
try
    [n, ts, sv] = plx_event_ts(filename, 257);
catch me
    if(strcmp(me.identifier,'MATLAB:unassignedOutputs'))
        error([mfilename ':NoTracking'],'Video tracking data not found in PLX file.');
    else rethrow(me);
    end
end

assert(n > 0, [mfilename ':NoTracking'],'Video tracking data not found in PLX file.');

if(exist('plx_vt_interpret','file')==2)
    [~, ~, mode, coords] = plx_vt_interpret(ts, sv);
else
    [coords, mode] = interpretPLXCoords(ts, sv);
end

end
