function [angle,coorddiffs] = GetAngleFromCoords(Coords)
%function GetAngleFromCoords gets the continuous bearing on a set of
%coordinates.  Feed in x and y matrix and get out a set of angles from pi
%to -pi that you can plot over time, or across time by average



if length(Coords(:,1))>1
    % get angle
    coorddiffs=[0,0;diff(Coords(:,2:3))];
    
    %
    vectorratio=(coorddiffs(:,2)./coorddiffs(:,1));
    % this is the angle of the bearing
    angle=atan(vectorratio);
    % now i have to figure out what to do if x and or y are <=0:
    for i=1:length(angle)
        % if x is negative and y is positive subtract 180 or pi
        if coorddiffs(i,1)<0 && coorddiffs(i,2)>=0
            % e.g. pi - my angle (it is negative so +-)
            angle(i)=angle(i)+pi;
            % if both x an y are negative, add 180 or pi
        elseif coorddiffs(i,1)<0 && coorddiffs(i,2)<=0
            angle(i)=angle(i)-pi;
            % if x is zero, you get nan and have to fix
        elseif coorddiffs(i,1)==0
            if(coorddiffs(i,2)>0)
                angle(i)=pi/2;
            elseif coorddiffs(i,2)<0
                angle(i)=-pi/2;
            end
        end
    end
else
    angle=nan; coorddiffs=nan;
end

