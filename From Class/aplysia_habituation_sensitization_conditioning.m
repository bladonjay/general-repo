%Simulate habituation in aplysia.
%
%INPUTS.
%  A = The decay
%  B 
%  x_s = The activity of the sensory neuron (mantle shelf).
%
%OUTPUT.
%  z_sm = The synaptic weight.

function [z_sm] = aplysia_habituation_sensitization_conditioning(A,x_s,B,x_f,C)

  N = length(x_s);      %The total number of time points.

  z_sm = zeros(N,1);      %Declare the output variables.
  
  z_sm(1) = 1;            %Set the initial conditions.
  
  dt = 0.01;
  for i=1:N-1
      z_sm(i+1) = z_sm(i) + dt*(-A*x_s(i) + B*x_f(i) + C*x_s(i)*x_f(i));
      
  end

end
  
  