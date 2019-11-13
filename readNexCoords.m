function [coords, mode] = readNexCoords(nex, varargin)
% pulls in coordinate data from a nex file, edited from unannotated Ben Kraus?
% JH Bladon edit.
% INPUTS:
%   nex: nexfile, path, or struct
%OUTPUTS:
%   coords: an m timestamps vector cols 1=ts, 2,3,4,5=coordinates of LED's

% JHB 1/21/2016
coords=[];


ip = inputParser;
ip.KeepUnmatched = true;
ip.addOptional('verbose',true);
ip.addParameter('plexonfcn',true);




ip.parse(varargin{:});
%%%%% patch by jhb
verbose=ip.Results.verbose;
plexonfcn=ip.Results.plexonfcn;
%%%%%%%


% really nice way to pull out the variables from the IP struct
for i=fields(ip.Results)', eval([i{1} '=ip.Results.' i{1} ';']); end;


% If we were not passed a structure (such as that returned by readNexFile),
% then assume it is either an FID or a filename.
if(~isstruct(nex));
    nex = readNexFileM(nex);
end

% Find the position data stored in the NEX file. It is assuming that the
% data came from a Plexon system, which stores the position data as the
% strobed event number 257 or as NEX or PLX coordinates.
markernum = [];
strobedind = nan;
idx =  1;

mode=[];
markers = [];

% find the names of the markers that you might want to use

if(isfield(nex,'markers'))
    for ii = 1:length(nex.markers)
        if any(strfind(nex.markers{ii}.name,'AVI'))
            % grab any avi coords, say number of coords and save their place
            markers = [markers,{['AVI ' num2str(length(nex.markers{ii}.timestamps))]}];
            markernum = [markernum;ii];
        end
        if any(strfind(nex.markers{ii}.name,'DVT'))
            markers = [markers,{['DVT ' num2str(length(nex.markers{ii}.timestamps))]}];
            markernum = [markernum;ii];
        end
        if any(strfind(lower(nex.markers{ii}.name),'strobed'))
            strobedind = ii;
            markers = [markers,{'Strobed'}];
            markernum = [markernum;ii];
        end
        if any(strfind(nex.markers{ii}.name,'PLX'))
            markers = [markers,{['PLX ' num2str(length(nex.markers{ii}.timestamps))]}];
            markernum = [markernum;ii];
        end
    end
    
    checked = checkBox(markers);
    
    % if no markers, dont even
    if any(checked)
        markernum = markernum(checked);
    else
        % if you didnt choose any markers, dont get coords!
        return
    end
    
    
end

% save space for where well get our frame markers
framermarkerid = nan;

% if there are events in the file, pull frame markers
if(isfield(nex,'events'))
    for ii = 1:length(nex.events)
        if (any(strfind(nex.events{ii}.name,'Frame Marker'))) || ...
                (any(strfind(nex.events{ii}.name,'Strobe')))
            
            framermarkerid = ii;
        end
    end
end


% If you found the tracking data in markers, translate it from strings to
% numbers, or use the frame marker to start your matrix
if(~isnan(strobedind))
    % (if we want it to talk to us)
    if(verbose); fprintf(1,'Position data found in NEX file (%d).\n',strobedind); end
    
    % if we have strobed data, pull the ts from that
    ts = nex.markers{strobedind}.timestamps;
    
    if(verbose); fprintf(1,'Converting marker strings into position values.\n'); end
    
    if(iscellstr(nex.markers{strobedind}.values{1}.strings))
        strs = char(nex.markers{strobedind}.values{1}.strings)';
    elseif(ischar(nex.markers{strobedind}.values{1}.strings))
        strs = nex.markers{strobedind}.values{1}.strings';
    end
    strs(end+1,:) = 0;
    sv = sscanf(strs,'%d');
    if(exist('plx_vt_interpret','file')==2 && plexonfcn)
        [~, ~, mode, coords] = plx_vt_interpret(ts, sv);
    else
        [coords, mode] = interpretPLXCoords(ts, sv);
    end
elseif ~isnan(framermarkerid)
    % otherwise, pull the indices from our event framemarker
    coords(:,1) = nex.events{framermarkerid}.timestamps;
else
    error('readNexCoords:NoTracking','Must have timestamps in file');
    
end
% preallocate the room for x and y coordinates with nan's


% now fill in frames where there are frames to fill
if (any(cellfun(@any,regexp(markers(checked),'AVI'))) ||...
        any(cellfun(@any,regexp(markers(checked),'DVT'))) || ...
        any(cellfun(@any,regexp(markers(checked),'PLX'))))
    
    if(~any(isnan(markernum)))
        % this should snap all coordinates to existing frame markers and
        % allow for empty space to be let in.
        
        tsAll=coords(:,1);
        
        % if you only selected a single tracking marker
        if length(markernum) == 1
            coords(:,2:3) = nan;
            X1=cellfun(@(a) str2num(a),nex.markers{markernum(1)}.values{1}.strings);
            Y1=cellfun(@(a) str2num(a),nex.markers{markernum(1)}.values{2}.strings);
            [~,bin1]=histc(nex.markers{markernum(1)}.timestamps,tsAll);
            X1(bin1==0) = [];
            Y1(bin1==0) = [];
            bin1(bin1==0) = [];
            coords(bin1,2) = X1;
            coords(bin1,3) = Y1;
            
            % if you selected two tracking markers
        elseif length(markernum) ==2
            coords(:,2:5) = nan;
            
            X1=cellfun(@(a) str2num(a),nex.markers{markernum(1)}.values{1}.strings);
            Y1=cellfun(@(a) str2num(a),nex.markers{markernum(1)}.values{2}.strings);
            [~,bin1]=histc(nex.markers{markernum(1)}.timestamps,tsAll);
            
            
            X1(bin1==0) = [];
            Y1(bin1==0) = [];
            bin1(bin1==0) = [];
            
            X2=cellfun(@(a) str2num(a),nex.markers{markernum(2)}.values{1}.strings);
            Y2=cellfun(@(a) str2num(a),nex.markers{markernum(2)}.values{2}.strings);
            [~,bin2]=histc(nex.markers{markernum(2)}.timestamps,tsAll);
            
            X2(bin2==0) = [];
            Y2(bin2==0) = [];
            bin2(bin2==0) = [];
            
            coords(bin1,2) = X1;
            coords(bin1,3) = Y1;
            coords(bin2,4) = X2;
            coords(bin2,5) = Y2;
        else
            fprintf('you selected %d markers,select only one or two tracking markers', length(markernum));
        end
        
        
        
        
    else
        warning('readNexCoords:NoTracking','Must have both Strobed and AVI or DVT coords');
        return
    end
    
    
    
end

end