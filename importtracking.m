function [xydata]=importtracking(input,validpolygon,smoothstamps)
%IMPORTTRACKING Import and smooth tracking data
%   Accepts the following input: PLX, NEX or raw LED tracking
%
%validpolygon -- 2D array of x-y positions defining a polygon.  All points
%               within polygon are included.  Also accepts "true", which
%               enables a GUI for outlining the valid polygon.
%smoothstamps -- scalar. Number of timestamps that are smoothed. Each
%               timestamp is about .03 seconds (default=15, or almost half
%               a second

%Handle input
if nargin < 2
    validpolygon=false;
end
if nargin < 3
    smoothstamps=15;
end

%Test for input type
if ischar(input)
    %Test whether input is a file
    if exist(input, 'file')==2
        [~,~,datatype]=fileparts(input);
        switch lower(datatype)
            case '.plx'
                xydata=readPLXCoords(input);
            case '.nex'
                [nexStruct]=readNexFileM(input);
                LED=cell(numel(nexStruct.markers),1);
                allts=[];
                for m=1:numel(nexStruct.markers)
                    ts=nexStruct.markers{m, 1}.timestamps;
                    LED{m}=[ts cellfun(@str2num,nexStruct.markers{m, 1}.values{1, 1}.strings)...
                        cellfun(@str2num,nexStruct.markers{m, 1}.values{2, 1}.strings)];
                    allts=[allts; ts]; %#ok<AGROW>
                end
                allts=unique(allts);
                
                xydata=allts(:);
                for m=1:numel(LED)
                    xydata(ismember(allts,LED{m}(:,1)),[2*m 2*m+1])=LED{m}(:,[2 3]);
                end
                %(inserts zeros for timestamps where no LED was detected)
            otherwise
                warning('Input is not a valid file type (PLX or NEX).');
                return
        end
    else
        warning('Input is not a valid file.');
        return
    end
else
    %Test whether xydata is a 3 or 5 column tracking vector with all
    %increasing timestamps
    if isnumeric(input) && any(size(input,2)==[3 5]) && all(diff(input(:,1),1,1)>0)
        xydata=input;
    else
        warning('Input is not valid.  It is neither a string or valid tracking matrix');
        return
    end
end

if size(xydata,2)==3
    xydata(:,[2 3]); %make xydata 5 column for tracking clean up
end

%Remove tracking that falls outside the predefined "validpolygon"
if all(validpolygon(:))
    %Create polygon
    if ~isnumeric(validpolygon)
        %Get background image
        [filename,directory]=uigetfile('*.avi','Select the AVI for image to map tracking.'); %Get video file
        vid=VideoReader([directory filename]); %Create video object
        im=read(vid,10); %Read 10th frame
        
        %Use GUI to define polygon
        [polyX,polyY]=outlinemaze2(xydata(:,[2 3]),im);
        validpolygon=[polyX polyY];
    end
    
    %Determine which coordinates are valid
    xdata=xydata(:,2:2:size(xydata,2));
    ydata=xydata(:,3:2:size(xydata,2));
    validpos=inpolygon(xdata,ydata,validpolygon(:,1),validpolygon(:,2));
    
    %Replace invalid tracking w/ nans (each LED separately)
    xydata(~validpos,2:2:size(xydata,2))=nan; %xdata
    xydata(~validpos,3:2:size(xydata,2))=nan; %ydata
    %Remove tracking entries where all tracking is invalid
    xydata(all(~validpos,2),:)=[];
end

%Clean and smooth tracking data
xydata=clean_plx_coords2(xydata,smoothstamps);%(homogenizes the timestamps!)
end