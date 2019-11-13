function [ V, taxis, n_spikes,tau ] = LIF( Iinput, R, C, Vrest, Vth, Vreset)
%Simulate an LIF neuron. This is is accompanied by leaky_integrate_fire.
% This 'neuron' will run for 400 iterations.
%
%INPUTS:
%  Iinput = input current, constant current in nanoamps
%  R = resistance. (in Mohms, start with 10)
%  C = capacitance. (in nanoFarads, start with like .5)
%  Vrest = resting voltage. (generally -70)
%  Vth = threshold voltage. (-60 or -55)
%  Vreset = voltage reset value. (reset to rest, -70)
%
%OUTPUT:
%  V = voltage of IF model.
%  taxis = time axis of simulation.
%  n_spikes = # spikes detected.

dt=0.1;                       %Set the timestep.
V(1)=Vrest;                   %Set the initial condition.

n_spikes = 0;                 %The # of spikes detected.

Vstar = Iinput*R + Vrest;     %Define Vstar = target voltage.
tau = R*C;                    %Define tau = time constant


for k=1:400-1                                  %For each time step,
    V(k+1) = V(k) + dt*(-(V(k)-Vstar)/tau);     %...evaluate V at the next time step
    if V(k+1) > Vth                             %When we cross threshold
        V(k+1) = Vreset;                        %...set V to reset
        n_spikes = n_spikes+1;                  %...and add 1 to spike counter.
    end
end

taxis = (1:length(V))*dt;         %Return the time axis too.

end

