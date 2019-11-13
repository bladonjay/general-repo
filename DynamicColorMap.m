function ptColors = DynamicColorMap(ptsX,ptsY,normX,normY,indexIntoNorm,radiusLimit)
%Index into norm has to fit normX(indexIntoNorm) = ptsX
%Can leave normX and normY empty to not normlize by another occupancy vector

% index to norm is an index of same size ptsX, and radiusLimit is the
% radius size to count close dots.

%Make them rows for easier
ptsX = ptsX(:)';
ptsY = ptsY(:)';

%Get distances: arrangement is row is anchor point, Y is to this pt
for ptI = 1:length(ptsX)
    distances(ptI,:) = cell2mat(arrayfun(@(x,y) hypot(ptsX(ptI)-x,ptsY(ptI)-y),ptsX,ptsY,'UniformOutput',false));
end

maxDist = max(max(distances));
minDist = min(min(distances(distances>0)));
    
%radiusLimit = 1;
for ptI = 1:length(ptsX)
    ptsClose(ptI,1) = sum(distances(ptI,:) <= radiusLimit) - 1;
end

%Repeat for occupancy normalizing
if any(normX) && any(normY)
    %if any(indexIntoNorm)
    %    normX = normX(indexIntoNorm);
    %    normY = normY(indexIntoNorm);
    %end
    normX = normX(:)';
    normY = normY(:)';
    
    for ptJ = 1:length(normX)
        distancesNorm(ptJ,:) = cell2mat(arrayfun(@(x,y) hypot(normX(ptJ)-x,normY(ptJ)-y),normX,normY,'UniformOutput',false));
    end
    
    for ptJ = 1:length(normX)
        ptsCloseNorm(ptJ,1) = sum(distancesNorm(ptJ,:) <= radiusLimit) - 1;
    end
    
    if any(indexIntoNorm)
        ptsCloseNorm = ptsCloseNorm(indexIntoNorm);
    end
    
    %These should now be the same size
    ptsClose = ptsClose./ptsCloseNorm;
    
end

maxClose = max(ptsClose);
minClose = min(ptsClose);

hh = figure;
cc = colormap(jet);
close(hh);

boundaries = linspace(minClose-0.00001,maxClose-0.00001,64);

ptColors = zeros(length(ptsX),3);
for bdStops = 1:64
    thesePts = ptsClose > boundaries(bdStops);
    ptColors(thesePts,:) = repmat(cc(bdStops,:),sum(thesePts),1);
end
    
end
