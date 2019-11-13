function map = joystickmap()
    eventcodes = maze_events();
    
    map = [1,  eventcodes.DIG;...   % "X" button
           3,  eventcodes.NODIG;... % "B" button
           4,  eventcodes.TM;...    % "Y" button
           2,  eventcodes.TMS;...   % "A" button
           6,  eventcodes.SAMPLE;...% "R" button
           11,  eventcodes.DECSPEED;... 
           12,  eventcodes.INCSPEED;... 
           17,  eventcodes.DECSPEED;...
           16,  eventcodes.INCSPEED;...
           ];
end
