function [ idx, distance ] = findclosest2D ( x_points, y_points, qx, qy)
%[ idx, distance ] = findclosest2D ( x_points, y_points, qx, qy)
%x_points and y_points are the data you're looking for points in
%qx and qy are your query points, can be a vector
%idx is index in x_/y_points

idx = nan(1,length(qx));
distance = nan(1,length(qx));
for queryPoint = 1:length(qx)
    dist = [];
    for checkThis=1:length(x_points)
        dist(checkThis) =...
            sqrt( (qx(queryPoint)-x_points(checkThis))^2 ...
                + (qy(queryPoint)-y_points(checkThis))^2 );
    end
    [distance(queryPoint), idx(queryPoint)] = min(dist);
end    

end