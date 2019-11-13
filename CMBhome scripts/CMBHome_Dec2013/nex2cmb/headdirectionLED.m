function [hd, offset, LEDloc] = headdirectionLED(c)

% Assume there are exactly two LEDs
assert(size(c,2)>=5,[mfilename ':NotEnoughCoordinates'],...
    'Input must have coordinates for at least two LEDs.');
if(size(c,2)>5);
    warning([mfilename ':IgnoringCoordinates'],...
    'This function uses only the first two sets of coordiantes and ignores the rest.');
end

% Convert zero coordiantes into NaNs, so they aren't used elsewhere.
c((c(:,2)==0 & c(:,3)==0),2:3) = NaN;
c((c(:,4)==0 & c(:,5)==0),2:3) = NaN;

% Calculate dx and dy from the LED coordinates.
dx = c(:,2)-c(:,4);
dy = c(:,3)-c(:,5);

% Calculate the distance between the LEDs.
dxy = sqrt(dx.^2+dy.^2);

% Determine whether the LEDs are mounted on the headstage or on the
% microdrive cup. The LEDs mounted on the headstage are closer together.
freqgap = mode(dxy);
if(freqgap>20); LEDloc = 'cup'; bigthreshold = 50; smallthreshold = 15;
else LEDloc = 'headstage'; bigthreshold = 16; smallthreshold = 9;
end

% Remove locations where the LEDs are too far away from each other.
biggap = dxy>bigthreshold;
dx(biggap) = NaN;
dy(biggap) = NaN;
dxy(biggap) = NaN;

% Remove locations where the LEDs are too close to each other.
smallgap = dxy<smallthreshold;
dx(smallgap) = NaN;
dy(smallgap) = NaN;
% dxy(smallgap) = NaN;

% Convert from dx and dy to an angle.
hd = atan2(dy,dx);

%{
% Find the mean position:
x = nanmean(c(:,2:2:end),2);
y = nanmean(c(:,3:2:end),2);

% Guess the frame rate based on the timestamps
dt = mode(diff(c(:,1)));

% It is probably 30 fps, so lets check that
if(dt >= 0.033 && dt <= 0.034); dt = 1/30; end

% Estimate the velocity of the rat.
dxt = [0; diff(x)];
dyt = [0; diff(y)];
speed = sqrt(dxt.^2 + dyt.^2)./dt;
hdt = atan2(dyt,dxt);

% Select a range of non-small speeds
speedrange = quantile(speed,[.7 .8]);

% Find when the rat is moving relatively quickly, but not so quickly it is
% probably a mistake.
running = (speed>=speedrange(1) & speed<=speedrange(2));
running = running & ~isnan(hd);

% Find the angle between a line connecting the LEDs and the actual head direction.
f = fittype('asin(sin(a+o))','independent','a','coefficients','o');
ofit = fit(hd(running),hdt(running),f,'StartPoint',0);
offset = coeffvalues(ofit);
%}

switch(LEDloc)
    case 'cup'; offset = -pi/2;
    case 'headstage'; offset = -pi;
end

% Rotate the angle 90 degress clockwise because LEDs are perpendicular to
% the rats head.
hd = mod(hd+pi+offset,2*pi)-pi;

%{
figure;
line(hdt(~running)+pi,hd(~running)+pi,'Color',[.5 .5 .5],'LineStyle','none','Marker','.');
line(hdt( running)+pi,hd( running)+pi,'Color',[ 0  0  1],'LineStyle','none','Marker','.');
line([0 2*pi],[0 2*pi],'Color','red','LineWidth',3);
title(sprintf('Offset: %.1f',offset*180/pi));
%}

% Attach the timestamps onto the head direction.
hd = [c(:,1) hd];

end