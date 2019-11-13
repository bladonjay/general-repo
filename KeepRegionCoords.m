function [newcoords,deletedcoords] = KeepRegionCoords(coords,mazeimage)
% function [newcoords,deletedcoords] = KeepRegionCoords(coords,mazeimage)
%   Removes spikes that occur outside of the boxed region, bins data and
%   warns that your coords dont mach up with pixels.  you can have coords
%   that are fractions, but they have to correspond to the pixels correctly

% coords must be inside the image:
if max(coords(:,2))>size(mazeimage,2)
    warning('X Coords occur outside of image bounds');
end
if max(coords(:,3))>size(mazeimage,1)
    warning('Y Coords occur outside of image bounds');
end

% snap to whole numbers;
pixelcoords=round(coords(:,2:3));

% get the logcal of each coordinate pair
for i=1:size(coords(:,1))
    keepcoords(i)=mazeimage(pixelcoords(i,2),pixelcoords(i,3));
end
% and tabulate the results
newcoords=coords(keepcoords,:);

deletedcoords=newcoords(~keepcoords,:);

fprintf('%.2f percent of coords got removed \n',nanmean(keepcoords)*100);



end

