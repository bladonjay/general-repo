function [v, w, taxis] = FHN(v_start, w_start, I, a, b, tau )
%Simulate the FHN model.
%
%INPUTS:
%  v_start = the initial value for v.
%  w_start = the initial value for w.
%  I = input drive.
%  a = recovery parameter.
%  b = recovery parameter.
%  tau = time constant for recovery.
%
%OUTPUT:
%  v = membrane potential or "voltage" 
%  w = recovery variable
%  t = time axis (useful for plotting).

dt=0.1;                       %Set the timestep.

v(1)=v_start;                    %Set the initial condition for v.
w(1)=w_start;                   %Set the initial condition for w.

for t=1:5000-1                                           %For each time step,
    v(t+1) = v(t) + dt*( v(t) - 1/3*v(t)^3 - w(t) + I);  %... update voltage,
    w(t+1) = w(t) + dt*( (v(t) - a - b*w(t))/tau );      %... update recovery.
    
end

taxis = (1:length(v))*dt;     %Return the time axis too.

end

