%Simulate Levy and Desmond model.
%
%INPUTS.
%  A = Hebbian growth rate.
%  B = Post-synaptically gated decrease rate.
%  x_i = The activity of the pre-synaptic neuron.
%  x_j = The activity of the post-synaptic neuron.
%
%OUTPUT.
%  z_ij = The synaptic weight.

function [z_ij] = levy_desmond(A,B,x_i,x_j)

  N = length(x_i);      %The total number of time points.

  z_ij = zeros(N,1);      %Declare the output variables.
  
  z_ij(1) = 1;            %Set the initial conditions.
  
  dt = 0.01;
  for i=1:N-1
      z_ij(i+1) = z_ij(i) + dt*(A*x_i(i)*x_j(i) - B*x_j(i)*z_ij(i));
      
  end

end
  
  