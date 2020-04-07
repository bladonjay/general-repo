function [newcoords,segments] = EpochCoords(oldcoords,events)
% function [newcoords,segments] = EpochCoords(oldcoords,events)
% Input a vector or matrix of coord data, and a nx2 matrix of start and
% stop events

% JHB



% make sure we have two columns
if size(events,2)~=2
    error('Your event matrix needs two columns');
end

% initiate the vector of keep indices
keepts=zeros(size(oldcoords,1),1);

% build a vector of indices that includes those that are in our epochs
for i=1:size(events,1)
    % grab coords that are within our epoch
    incoords=oldcoords(:,1)>events(i,1) & oldcoords(:,1)<events(i,2);
    % contribute to talley
    keepts=keepts | incoords;
end

% newcoords are the rows of any ts thats within our epochs
newcoords=oldcoords(logical(keepts),:);
% and segments are when those transitions occur.
segments=find(diff(keepts)~=0);

end

