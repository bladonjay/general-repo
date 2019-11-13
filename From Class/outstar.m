function [x,y,z] = outstar(N,I,J, A,B,C,D,E,F,G)
% N= num timesteps
% I= number of cells
% J= Initial conditions
% A=Decay rate of sensory neuron activities.
% B=caling of conditioned stimulus input to sensory neurons.
% C=Decay rate of motor neuron activities.
% D=Scaling of input from sensory to motor neurons.
% E=Scaling of unconditioned stimulus input to motor neurons.
% F=Hebbian learning rate.
% G=Decay rate of synaptic weights.
  ncells_x=3;         %The # of sensory neurons.
  
  x = zeros(N,ncells_x); % sensory
  ncells_y=3;
  y = zeros(N,ncells_y); % motor
  z = zeros(N,ncells_x,ncells_y); % weights
  
  x(1,:) = zeros(1,ncells_x); %
  y(1,:) = zeros(1,ncells_y);
  z(1,:,:)=ones(ncells_x,ncells_y);
  
  dt = 0.1;               %Set the timestep.
  
  for t=1:N-1                             %For each time step.
      z0 = squeeze(z(t,:,:));             %... determine the synaptic weights now.
      
      for i=1:ncells_x                    %... step forward the sensory neuron activities.
          % A is sensory decay, B= sensory weight
          x(t+1,i) = x(t,i) + dt*(-A*x(t,i) + B*I(i));
      end
      
      for j=1:ncells_y                    %... step forward the motor neuron activities.
          % C= motor decay, D= sensory to motor cxn E=motor input scale
          y(t+1,j) = y(t,j) + dt*(-C*y(t,j) + D*sum(x(t,:).*z(t,:,j)) + E*J(j));
      end
      
      for i=1:ncells_x                    %... step forward the synaptic weight activities.
          for j=1:ncells_y
              % G= decay of all weights, f=hebbian learning rate
              % synaptic decay in outstar depends on sensory decay
              z(t+1,i,j) = z(t,i,j) + dt*(x(t,i)*(-G*z(t,i,j) + F*y(t,j)));
          end
      end
  end
  
end
  

  
  