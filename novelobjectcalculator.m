function [obj1tdist,obj2tdist,obj1sampletime,obj2sampletime,finalratio]=novelobjectcalculator(session,obj1coords,obj2coords,varargin)

% just a quick calculator to calc differential exploration time if oyu have
% coord data and the coordinates of the objects

obj1tdist=[]; obj2tdist=[]; obj1sampletime=[]; obj2sampletime=[]; finalratio=[];


ip=inputParser;
addOptional(ip,'sampledist',100);
addOptional(ip,'minutestokeep',3);
addOptional(ip,'plotter',0);
parse(ip,varargin{:});

sampledist=ip.Results.sampledist;
minutestokeep=ip.Results.minutestokeep;
plotter=ip.Results.plotter;

coords=session.edit_coords(1:minutestokeep*60*30,:);
coords(:,2)=smooth(coords(:,2)); coords(:,3)=smooth(coords(:,3));



% calculate how far the rat was from each oject and how long he was close
obj1time(:,1)=coords(:,2)-obj1coords(1);
obj1time(:,2)=coords(:,3)-obj1coords(2);
obj1tdist=sqrt(obj1time(:,1).^2+obj1time(:,2).^2);
obj1sampleframes=obj1tdist<sampledist;  % 3 pixels, or 15 cm
obj1sampletime=sum(obj1sampleframes*.033);

obj2time(:,1)=coords(:,2)-obj2coords(:,1);
obj2time(:,2)=coords(:,3)-obj2coords(:,2);
obj2tdist=sqrt(obj2time(:,1).^2+obj2time(:,2).^2);
obj2sampleframes=obj2tdist<sampledist;
obj2sampletime=sum(obj2sampleframes*.033);
if plotter
    figure;
    plot(coords(:,1)',obj1tdist,coords(:,1)',obj2tdist,'r');
    legend('Object 1 Distance','Object 2 Distance');
end

% object 1- object 2 over total object exploration
finalratio=(obj1sampletime-obj2sampletime) / (obj1sampletime+obj2sampletime);

end