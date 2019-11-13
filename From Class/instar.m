function [x,y,z] = instar(N,I,J, A,B,C,D,E,F,G)
% 
%
% INPUTS
%   N= number of timepoints
%   I= input activity to start
%   J= output activity to start

  ncells_x=length(I);         %The # of sensory neurons.
  x = zeros(N,ncells_x);
  ncells_y=length(J);
  y = zeros(N,ncells_y);
  z = zeros(N,ncells_x,ncells_y);
  
  % initialize cells to start off
  x(1,:) = zeros(1,ncells_x);
  y(1,:) = zeros(1,ncells_y);
  % initalize synaptic weights to start at 1
  z(1,:,:)=ones(ncells_x,ncells_y);
  
  dt = 0.1;               %Set the timestep.
  
  for t=1:N-1                             %For each time step.
      %z0 = squeeze(z(t,:,:));             %... determine the synaptic weights now.
      
      for i=1:ncells_x                    %... step forward the sensory neuron activities.
          x(t+1,i) = x(t,i) + dt*(-A*x(t,i) + B*I(i));
      end
      
      for j=1:ncells_y                    %... step forward the motor neuron activities.
          y(t+1,j) = y(t,j) + dt*(-C*y(t,j) + D*sum(x(t,:).*z(t,:,j)) + E*J(j));
      end
      
      for i=1:ncells_x                    %... step forward the synaptic weight activities.
          for j=1:ncells_y
              % synaptic decay in instar depends on motor activity
              % G= decay, F= hebbian learning
              z(t+1,i,j) = z(t,i,j) + dt*(y(t,j)*(-G*z(t,i,j) + F*x(t,i)));
          end
      end
  end
  
end
  

  
  