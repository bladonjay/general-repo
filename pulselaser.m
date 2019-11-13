function pulselaser(sport, delay, pulse)
% delay and pulse are duration in seconds
% maximum value for each are 2^16-1 = 65535ms = 65.535s
% but thats only because we pulse 2 bytes and our numbers are in ms
% the baudrate here is 4khz so you could bring this down to .25 ms if you
% were fast enough (maybe good for stim, then you could do duty cycle)
% 

delay = delay*1000;
pulse = pulse*1000;

assert(isa(sport,'serial'),'pulselaser:invalidserial','First input argument must be a handle to a serial port.');
assert(delay>=0,  'pulselaser:negativedelay','Delay length must be greater than or equal to zero');
assert(pulse>=0,  'pulselaser:negativedelay','Delay length must be greater than or equal to zero');
assert(delay<2^16,'pulselaser:delaytoolong','Delay length must be less than 65.535s');
assert(pulse<2^16,'pulselaser:pulsetoolong','Delay length must be less than 65.535s');
assert(strcmp(sport.Status,'open'),'pulselaser:portclosed','Serial port must be opened before calling pulselaser');

% He must have programmed the arduino to listen to the serial port and
% pulse the laser based on the two bytes

dhbyte = uint8(floor(delay/256)); %
dlbyte = uint8(mod(delay,256));
phbyte = uint8(floor(pulse/256));
plbyte = uint8(mod(pulse,256));

fwrite(sport, [dhbyte, dlbyte, phbyte, plbyte]);

end