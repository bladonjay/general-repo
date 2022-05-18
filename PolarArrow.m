function [hfig] = PolarArrow(rho,resultant,hfig,mycolor)
% function [hfig] = PolarArrow(rho,resultant,hfig,mycolor)
% generates arrow in current polarplot;


if ~exist('hfig','var')
    hfig=gcf;
end
if isempty(hfig)
    hfig=gcf;
end

if ~exist('mycolor','var')
    mycolor=rgbcolormap('red');
end


%rho = rand(1)*2*pi;
%resultant = 0.5 + (1-0.5).*rand(1);
%%%%arrow head %%%%
arrowhead_length    = resultant/5; % arrow head length relative to resultant_length
num_arrowlines = 100;
arrowhead_angle = deg2rad(30); % degrees
%%%%arrow tip coordinates %%%%
t1 = repmat(rho,1,num_arrowlines);
r1 = repmat(resultant,1,num_arrowlines);
%%%%arrow base coordinates %%%%
b = arrowhead_length.*tan(linspace(0,arrowhead_angle,num_arrowlines/2));
theta = atan(b./(resultant-arrowhead_length));
pre_t2 = [theta, -theta];
r2 = (resultant-arrowhead_length)./cos(pre_t2);
t2 = t1(1)+pre_t2;
%%%%plot %%%%
polarplot([t1(1) t1(1)],[0 r1(1)-0.9*arrowhead_length],'color',mycolor,'linewidth',3)
hold on
polarplot([t1; t2],[r1; r2],'color',mycolor);
end

