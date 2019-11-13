function [smoothed, interprows] = CleanCoordData(coorddata,smoothtype,gap)
%  function [smoothdata, interprows] = CleanCoordData(coorddata,smoothtype,ga[)
% interpolates and smooths coordinate data that has nans, or will just
% smooth if no nans
% INPUTS:
% coorddata: nx3 matrix, [ts coord1 coord2]
% smoothtype: ether a span (For linear) or a method
% gap: maximum gap between timestamps to interpolate & smooth across

cutjumps=1;

smoothed=[]; interprows=[];
%% figure out how we're going to smooth our data
if ~exist('smoothtype','var')
    smoothtype='movmean';
    span=7;
else
    if isnumeric(smoothtype)
        span=smoothtype;
        smoothtype='movmean';
    end
end
%% make sure our gaps are big enough
if ~exist('gap','var')
    gap=1;
end

% squish into one xy coord
if size(coorddata,2)>3
    coorddata(:,3)=nanmean(coorddata(:,[3 5]),2);
    coorddata(:,2)=nanmean(coorddata(:,[2 4]),2);
    coorddata=coorddata(:,1:3);
end

%% now find out where our contiguous chunks of data are
% smooths each chunk of data where the ts are less than 1 second apart
% the find will index the last timestamp of each block
chunks=[0 find(diff(coorddata(:,1)>gap)) length(coorddata(:,1))];

%% gotta go through and erase when the coords skip super far

firstcoord=find(~isnan(coorddata(:,2)),1,'first');
lastcoord=coorddata(firstcoord,:);
%
if cutjumps
displacement=hypot(diff(coorddata(:,2)),diff(coorddata(:,3)));
jumpthresh=prctile(displacement./diff(coorddata(:,1)),99);

%
for i=firstcoord+1:length(coorddata)
    % get the distance between that and the next timestamp
    coordjump=hypot(lastcoord(2)-coorddata(i,2),lastcoord(3)-coorddata(i,3))...
        /(coorddata(i,1)-lastcoord(1));
    timejump=coorddata(i,1)-lastcoord(1);
    if ~isnan(coordjump) && timejump<.1
        % if its real and its close in time
        if coordjump>jumpthresh
            coorddata(i,[2 3])=[nan nan];
        end
    end
    % if that coord survived, its the new reference spike
    if ~isnan(coorddata(i,2))
        lastcoord=coorddata(i,:);
    end
end
end

%% for each chunk, interpolate and smooth
for i=1:length(chunks)-1
    
    mychunk=coorddata(chunks(i)+1:chunks(i+1),:);
    nanrows=isnan(mychunk(:,2));
    % if any nans, interpolate them
    try
        if sum(nanrows)>0
            mychunk(nanrows,2:end)=interp1(mychunk(~nanrows,1),mychunk(~nanrows,2:end),mychunk(nanrows,1));
        end
        smoothchunk=mychunk(:,1);
        % and smooth a bit
        smoothchunk(:,2)=smoothdata(mychunk(:,2),smoothtype,span);
        smoothchunk(:,3)=smoothdata(mychunk(:,3),smoothtype,span);
        
        % now tack on to your full dataset
        smoothed=[smoothed; smoothchunk];
        interprows=[interprows; nanrows];
    catch
        warning('Chunk %.2f to %.2f not included in Coord Cleanup \n',chunks(i),chunks(i+1));
    end
end

end
