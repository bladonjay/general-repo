function [coords, mode] = readNexCoords(nex, varargin)

ip = inputParser;
ip.KeepUnmatched = true;
ip.addParamValue('verbose',false);
ip.addParamValue('plexonfcn',true);



ip.parse(varargin{:});

for i=fields(ip.Results)', eval([i{1} '=ip.Results.' i{1} ';']); end;


% If we were not passed a structure (such as that returned by readNexFile),
% then assume it is either an FID or a filename.
if(~isstruct(nex)); 
    nex = readNexFileM(nex);
end

% Find the position data stored in the NEX file. It is assuming that the
% data came from a Plexon system, which stores the position data as the
% strobbed event number 257 or as NEX or PLX coordinates.
trackingind = [];
strobedind = nan;
idx =  1;


markers = [];

% find the names of the markers that you might want to use

if(isfield(nex,'markers'))
 
    for ii = 1:length(nex.markers)
        
        
        if any(strfind(nex.markers{ii}.name,'AVI'))
           
            markers = [markers,{['AVI ' num2str(length(nex.markers{ii}.timestamps))]}];
            trackingind = [trackingind;ii];
        end
        
        if any(strfind(nex.markers{ii}.name,'DVT'))
           
            markers = [markers,{['DVT ' num2str(length(nex.markers{ii}.timestamps))]}];
             trackingind = [trackingind;ii];
        end
        
        if any(strfind(nex.markers{ii}.name,'Strobed'))
          strobedind = ii;
            markers = [markers,{'Strobed'}];
             trackingind = [trackingind;ii];
        end
        
        if any(strfind(nex.markers{ii}.name,'PLX'))
          
             markers = [markers,{['PLX ' num2str(length(nex.markers{ii}.timestamps))]}];
             trackingind = [trackingind;ii];
        end
        
    end
       checked = checkBox(markers);
trackingind = trackingind(checked);


end
framermarkerid = nan;

% if there are events in the file, pull frame markers
if(isfield(nex,'events'))
    for ii = 1:length(nex.events)
        
        if (any(strfind(nex.events{ii}.name,'Frame Marker')))
            framermarkerid = ii;
            
        end
        
    end
end


% If you found the tracking data in events, translate it from strings to
% numbers, or use the frame marker to start your matrix
if(~isnan(strobedind))
    if(verbose); fprintf(1,'Position data found in NEX file (%d).\n',strobedind); end
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
    coords(:,1) = nex.events{framermarkerid}.timestamps;
else
    error('readNexCoords:NoTracking','Must have timestamps in file');
    
end
% preallocate the room for x and y coordinates with nan's


% now fill in frames where there are frames to fill
if (any(cellfun(@any,regexp(markers(checked),'AVI'))) ||...
        any(cellfun(@any,regexp(markers(checked),'DVT'))) || ...
        any(cellfun(@any,regexp(markers(checked),'PLX'))));
    
    if(~any(isnan(trackingind)))
         % this should snap all coordinates to existing frame markers and
        % allow for empty space to be let in.
        
        tsAll=coords(:,1);
        
        % if you only selected a single tracking marker
        if length(trackingind) == 1 
            coords(:,2:3) = nan;
            X1=cellfun(@(a) str2num(a),nex.markers{trackingind(1)}.values{1}.strings);
            Y1=cellfun(@(a) str2num(a),nex.markers{trackingind(1)}.values{2}.strings);
            [~,bin1]=histc(nex.markers{trackingind(1)}.timestamps,tsAll);
            X1(bin1==0) = [];
            Y1(bin1==0) = [];
            bin1(bin1==0) = [];
            coords(bin1,2) = X1;
            coords(bin1,3) = Y1;
            
        % if you selected two tracking markers   
        elseif length(trackingind) ==2
            coords(:,2:5) = nan;
            
            X1=cellfun(@(a) str2num(a),nex.markers{trackingind(1)}.values{1}.strings);
            Y1=cellfun(@(a) str2num(a),nex.markers{trackingind(1)}.values{2}.strings);
            [~,bin1]=histc(nex.markers{trackingind(1)}.timestamps,tsAll);
            
            
            X1(bin1==0) = [];
            Y1(bin1==0) = [];
            bin1(bin1==0) = [];
            
            X2=cellfun(@(a) str2num(a),nex.markers{trackingind(2)}.values{1}.strings);
            Y2=cellfun(@(a) str2num(a),nex.markers{trackingind(2)}.values{2}.strings);
            [~,bin2]=histc(nex.markers{trackingind(2)}.timestamps,tsAll);
            
            X2(bin2==0) = [];
            Y2(bin2==0) = [];
            bin2(bin2==0) = [];
            
            coords(bin1,2) = X1;
            coords(bin1,3) = Y1;
            coords(bin2,4) = X2;
            coords(bin2,5) = Y2;
        else
            fprintf('you selected %d markers,select only one or two tracking markers', trackingind);
        end
        
        
        
        
    else
        error('readNexCoords:NoTracking','Must have both Strobed and AVI coords');
    end
    
    
    
end

end