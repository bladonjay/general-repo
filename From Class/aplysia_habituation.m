%Simulate habituation in aplysia.
%
%INPUTS.
%  A = The strength of habituation.
%  x_s = The activity of the sensory neuron (mantle shelf).
%
%OUTPUT.
%  z_sm = The synaptic weight.

function [z_sm] = aplysia_habituation(A,x_s)

  N = length(x_s);        %The total number of time points.

  z_sm = zeros(N,1);      %Declare the output variables.
  
  z_sm(1) = 1;            %Set the initial conditions.
  
  % dt is like 10 ms
  dt = 0.01;
  % for each step shrink the synapse if it fires
  for i=1:N-1
      z_sm(i+1) = z_sm(i) + dt*(-A*x_s(i));
      
  end

end
  
  