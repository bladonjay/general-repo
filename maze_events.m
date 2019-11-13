function [eventcodes] = maze_events()
    % These event codes are hard-wired inputs from the button-box. The left
    % side is what we have named our buttons on the machine, and the right
    % side is the number that we've added so that it can go into the map
    % box as a code
    % I think the server only accepts a few of these
    
    
%     eventcodes.LASERON = 3;      % Laser turned on
%     eventcodes.LASEROFF = 4;      % Laser turned off.
%     eventcodes.TL = 3;      % Top Left button
%     eventcodes.TR = 4;      % Top Right button
%     eventcodes.ML = 5;      % Middle Left button
%     eventcodes.MR = 6;      % Middle Right button
%     eventcodes.BL = 7;      % Bottom Left button
%     eventcodes.BR = 8;      % Bottom Right button
 
    % These event codes are hard-wired from other devices.
%     eventcodes.RVALVE = 9;  % The right valve opening
%     eventcodes.LVALVE = 10; % The left valve opening
%     eventcodes.DACLD  = 11; % DAC updated to new value
%     eventcodes.SPEED1 = 12; % First speed sensor
%     eventcodes.SPEED2 = 13; % Second speed sensor
%     eventcodes.MVALVE = 14; % The middle valve opening
%     eventcodes.EVALVE = 15; % The extra valve opening

    
    % These event codes just set a convention for recording eventcodes.
    eventcodes.LASERZONE  = 101;% Laser triggered
%     eventcodes.LL  = 102;	% Entered the left reward area using a valid path.
%     eventcodes.RW  = 103;	% Entered the right reward area using an invalid path.
%     eventcodes.LW  = 104;	% Entered the left reward area using an invalid path.
    eventcodes.DIG  = 105;  % Rat dug in the pot.
    eventcodes.NODIG = 106; % Rat did not dig in the pot.
    eventcodes.TM  = 107;	% Entered the treadmill from the correct direction.
    eventcodes.TMS = 108;   % Button press used to stop the treadmill.
    eventcodes.SAMPLE  = 109;   % Rat sampled object or/and odor.
    eventcodes.INCSPEED =110;
    eventcodes.DECSPEED =111;
    % added by JH Bladon
    eventcodes.BOARDSTART = 12;
    eventcodes.BOARDEND = 13;
    eventcodes.DOOR1=14;
    eventcodes.DOOR2=15;
    eventcodes.QuickTimer=800; % just need a placeholder here
    
    
    % MAP System (hard coded)
    eventcodes.CINEPLEX = 257;
    eventcodes.START = 258;
    eventcodes.STOP = 259;
end
