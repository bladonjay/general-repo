%% Assignment 11 edits


%% Q1:4
N=1000;             %The total # of time steps.

ncells_x=3;         %The # of sensory neurons.
ncells_y=3;

I = ones(1,ncells_x);
J = zeros(1,ncells_y);  
J(:,1)=0.1;
J(:,2)=0.4;
J(:,3)=0.9;

A=1;        %Decay rate of sensory neuron activities.
B=1;        %Scaling of conditioned stimulus input to sensory neurons.
C=1;        %Decay rate of motor neuron activities.
D=1;        %Scaling of input from sensory to motor neurons.
E=0.5;      %Scaling of unconditioned stimulus input to motor neurons.
F=0.1;      %Hebbian learning rate.
G=1;        %Decay rate of synaptic weights.

[x,y,z] = instar(N,I,J, A,B,C,D,E,F,G);
figure(1)
plot(x)
xlabel('Time [steps]');  ylabel('Sensory neuron activities (x)')
legend({'x_1'; 'x_2'; 'x_3'})

figure(2)
plot(y)
xlabel('Time [steps]');  ylabel('Motor neuron activities (y)')
legend({'y_1'; 'y_2'; 'y_3'})

figure(3);                  %Make an empty figure,
color = {'r'; 'g'; 'b'};    %Define colors for line plots.
shape = {'-'; '--'; ':'};   %Define symbols for line plots.
clf();                      %Clear the figure.
for i=1:ncells_x            %For each sensory neuron,
    for j=1:ncells_y        %... and for each motor neuron,
        hold on             %... hold on to the graphics window,
        plot(z(:,i,j), [color{j} shape{i}])  %... and plot the synaptic weight dynamics, with a chosen color and shape.
        hold off            %... then release the graphics window,
    end                     %... and close the for-loops.
end
xlabel('Time [steps]');  ylabel('Synaptic weights (z)')
legend({'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'});


% sensory activities are all 1, motor actiities and symaptic weights all go
% past 10
% this is becaues the synaptic plasticity dominates because everything is 1

%% Q5: 7
N=1000;             %The total # of time steps.

ncells_x=3;         %The # of sensory neurons.
ncells_y=3;

I = ones(1,ncells_x);

J = zeros(1,ncells_y);  
J(:,1)=0.1;
J(:,2)=0.4;
J(:,3)=0.9;

A=1;        %Decay rate of sensory neuron activities.
B=1;        %Scaling of conditioned stimulus input to sensory neurons.
C=1;        %Decay rate of motor neuron activities.
D=1;        %Scaling of input from sensory to motor neurons.
E=0.5;      %Scaling of unconditioned stimulus input to motor neurons.
F=0.1;      %Hebbian learning rate.
G=1;        %Decay rate of synaptic weights.

[x,y,z] = outstar(N,I,J, A,B,C,D,E,F,G);
figure(1)
plot(x)
xlabel('Time [steps]');  ylabel('Sensory neuron activities (x)')
legend({'x_1'; 'x_2'; 'x_3'})

figure(2)
plot(y)
xlabel('Time [steps]');  ylabel('Motor neuron activities (y)')
legend({'y_1'; 'y_2'; 'y_3'})

figure(3);                  %Make an empty figure,
color = {'r'; 'g'; 'b'};    %Define colors for line plots.
shape = {'-'; '--'; ':'};   %Define symbols for line plots.
clf();                      %Clear the figure.
for i=1:ncells_x            %For each sensory neuron,
    for j=1:ncells_y        %... and for each motor neuron,
        hold on             %... hold on to the graphics window,
        plot(z(:,i,j), [color{j} shape{i}])  %... and plot the synaptic weight dynamics, with a chosen color and shape.
        hold off            %... then release the graphics window,
    end                     %... and close the for-loops.
end
xlabel('Time [steps]');  ylabel('Synaptic weights (z)')
legend({'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'});

% sensory all go to 1, motor all go between 0 and 1, and synaptic weights
% all go 0.0071, .0286, 0.0643
%% 
% q8 the decay term that determines the synaptic rates is multiplied by the
% sensory activity

% q9 you swap the x for the y now to determine the synaptic decay term
%% Q10

N=1000;             %The total # of time steps.

ncells_x=3;         %The # of sensory neurons.
ncells_y=3;

I = ones(1,ncells_x);

J = zeros(1,ncells_y);
J(:,1)=1;
J(:,2)=0;
J(:,3)=0;

A=1;        %Decay rate of sensory neuron activities.
B=1;        %Scaling of conditioned stimulus input to sensory neurons.
C=1;        %Decay rate of motor neuron activities.
D=1;        %Scaling of input from sensory to motor neurons.
E=0.5;      %Scaling of unconditioned stimulus input to motor neurons.
F=0.1;      %Hebbian learning rate.
G=1;        %Decay rate of synaptic weights.

[x,y,z] = instar(N,I,J, A,B,C,D,E,F,G);
figure(1)
plot(x)
xlabel('Time [steps]');  ylabel('Sensory neuron activities (x)')
legend({'x_1'; 'x_2'; 'x_3'})

figure(2)
plot(y)
xlabel('Time [steps]');  ylabel('Motor neuron activities (y)')
legend({'y_1'; 'y_2'; 'y_3'})

figure(3);                  %Make an empty figure,
color = {'r'; 'g'; 'b'};    %Define colors for line plots.
shape = {'-'; '--'; ':'};   %Define symbols for line plots.
clf();                      %Clear the figure.
for i=1:ncells_x            %For each sensory neuron,
    for j=1:ncells_y        %... and for each motor neuron,
        hold on             %... hold on to the graphics window,
        plot(z(:,i,j), [color{j} shape{i}])  %... and plot the synaptic weight dynamics, with a chosen color and shape.
        hold off            %... then release the graphics window,
    end                     %... and close the for-loops.
end
xlabel('Time [steps]');  ylabel('Synaptic weights (z)')
legend({'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'});

% the sensory activities converge to 1, the motor activities are between 0
% and 1 (0.8 and about 0.3) and the synaptic weights end up at 0.1 for all

%% Q 14
N=1000;             %The total # of time steps.

ncells_x=3;         %The # of sensory neurons.
ncells_y=3;

I = zeros(1,ncells_x); I(:,1)=1; I(:,2)=0; I(:,3)=0;  
J = zeros(1,ncells_y); J(:,1)=1; J(:,2)=0; J(:,3)=0; 


A=1;        %Decay rate of sensory neuron activities.
B=1;        %Scaling of conditioned stimulus input to sensory neurons.
C=1;        %Decay rate of motor neuron activities.
D=1;        %Scaling of input from sensory to motor neurons.
E=0.5;      %Scaling of unconditioned stimulus input to motor neurons.
F=0.1;      %Hebbian learning rate.
G=1;        %Decay rate of synaptic weights.

[x,y,z] = instar(N,I,J, A,B,C,D,E,F,G);
figure(1)
plot(x)
xlabel('Time [steps]');  ylabel('Sensory neuron activities (x)')
legend({'x_1'; 'x_2'; 'x_3'})

figure(2)
plot(y)
xlabel('Time [steps]');  ylabel('Motor neuron activities (y)')
legend({'y_1'; 'y_2'; 'y_3'})

figure(3);                  %Make an empty figure,
color = {'r'; 'g'; 'b'};    %Define colors for line plots.
shape = {'-'; '--'; ':'};   %Define symbols for line plots.
clf();                      %Clear the figure.
legendnames={'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'};
for i=1:ncells_x            %For each sensory neuron,
    for j=1:ncells_y        %... and for each motor neuron,
        subplot(3,3,3*(i-1)+j)
        hold on             %... hold on to the graphics window,
        plot(z(:,i,j), [color{j} shape{i}])  %... and plot the synaptic weight dynamics, with a chosen color and shape.
        legend(legendnames{3*(i-1)+j});
        hold off            %... then release the graphics window,
        
    end                     %... and close the for-loops.
end
xlabel('Time [steps]');  ylabel('Synaptic weights (z)')
legend({'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'});
% all weights from s1 end up at 1, whereas the others are at 0
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Assignment 12 %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Q 1
N=1000;             %The total # of time steps.

ncells_x=3;         %The # of sensory neurons.
ncells_y=3;

I = zeros(1,ncells_x); I(:,1)=1; I(:,2)=1; I(:,3)=0;  
J = zeros(1,ncells_y); J(:,1)=0.9; J(:,2)=0.6; J(:,3)=0.3; 


A=1;        %Decay rate of sensory neuron activities.
B=1;        %Scaling of conditioned stimulus input to sensory neurons.
C=1;        %Decay rate of motor neuron activities.
D=1;        %Scaling of input from sensory to motor neurons.
E=0.5;      %Scaling of unconditioned stimulus input to motor neurons.
F=0.1;      %Hebbian learning rate.
G=1;        %Decay rate of synaptic weights.

%%

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
%%
figure(1)
plot(x)
xlabel('Time [steps]');  ylabel('Sensory neuron activities (x)')
legend({'x_1'; 'x_2'; 'x_3'})

figure(2)
plot(y)
xlabel('Time [steps]');  ylabel('Motor neuron activities (y)')
legend({'y_1'; 'y_2'; 'y_3'})

figure(3);                  %Make an empty figure,
color = {'r'; 'g'; 'b'};    %Define colors for line plots.
shape = {'-'; '--'; ':'};   %Define symbols for line plots.
clf();                      %Clear the figure.
legendnames={'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'};
for i=1:ncells_x            %For each sensory neuron,
    for j=1:ncells_y        %... and for each motor neuron,
        subplot(3,3,3*(i-1)+j)
        hold on             %... hold on to the graphics window,
        plot(z(:,i,j), [color{j} shape{i}])  %... and plot the synaptic weight dynamics, with a chosen color and shape.
        legend(legendnames{3*(i-1)+j});
        hold off            %... then release the graphics window,
        
    end                     %... and close the for-loops.
end
xlabel('Time [steps]');  ylabel('Synaptic weights (z)')
%legend({'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'});

%% Q 4
N=1000;             %The total # of time steps.

ncells_x=3;         %The # of sensory neurons.
ncells_y=3;

I = zeros(1,ncells_x); I(:,1)=.1; I(:,2)=.4; I(:,3)=0.9;  
J = zeros(1,ncells_y); J(:,1)=.1; J(:,2)=0.4; J(:,3)=0.9; 


A=1;        %Decay rate of sensory neuron activities.
B=1;        %Scaling of conditioned stimulus input to sensory neurons.
C=1;        %Decay rate of motor neuron activities.
D=1;        %Scaling of input from sensory to motor neurons.
E=0.5;      %Scaling of unconditioned stimulus input to motor neurons.
F=0.1;      %Hebbian learning rate.
G=1;        %Decay rate of synaptic weights.

[x,y,z] = instar(N,I,J, A,B,C,D,E,F,G);
figure(1)
plot(x)
xlabel('Time [steps]');  ylabel('Sensory neuron activities (x)')
legend({'x_1'; 'x_2'; 'x_3'})

figure(2)
plot(y)
xlabel('Time [steps]');  ylabel('Motor neuron activities (y)')
legend({'y_1'; 'y_2'; 'y_3'})

figure(3);                  %Make an empty figure,
color = {'r'; 'g'; 'b'};    %Define colors for line plots.
shape = {'-'; '--'; ':'};   %Define symbols for line plots.
clf();                      %Clear the figure.
legendnames={'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'};
for i=1:ncells_x            %For each sensory neuron,
    for j=1:ncells_y        %... and for each motor neuron,
        subplot(3,3,3*(i-1)+j)
        hold on             %... hold on to the graphics window,
        plot(z(:,i,j), [color{j} shape{i}])  %... and plot the synaptic weight dynamics, with a chosen color and shape.
        legend(legendnames{3*(i-1)+j});
        hold off            %... then release the graphics window,
        
    end                     %... and close the for-loops.
end
xlabel('Time [steps]');  ylabel('Synaptic weights (z)')
%legend({'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'});

%% q 6
N=1000;             %The total # of time steps.

ncells_x=3;         %The # of sensory neurons.
ncells_y=3;

I = zeros(1,ncells_x); 
I(:,1)=1;
I(:,2)=1;
I(:,3)=0;
J = zeros(1,ncells_y);
J(:,1)=0;
J(:,2)=0;
J(:,3)=1;



A=1;        %Decay rate of sensory neuron activities.
B=1;        %Scaling of conditioned stimulus input to sensory neurons.
C=1;        %Decay rate of motor neuron activities.
D=1;        %Scaling of input from sensory to motor neurons.
E=0.1;      %Scaling of unconditioned stimulus input to motor neurons.
F=0.5;      %Hebbian learning rate.
G=1;        %Decay rate of synaptic weights.

[x,y,z] = instar(N,I,J, A,B,C,D,E,F,G);
figure(1)
plot(x)
xlabel('Time [steps]');  ylabel('Sensory neuron activities (x)')
legend({'x_1'; 'x_2'; 'x_3'})

figure(2)
plot(y)
xlabel('Time [steps]');  ylabel('Motor neuron activities (y)')
legend({'y_1'; 'y_2'; 'y_3'})

figure(3);                  %Make an empty figure,
color = {'r'; 'g'; 'b'};    %Define colors for line plots.
shape = {'-'; '--'; ':'};   %Define symbols for line plots.
clf();                      %Clear the figure.
legendnames={'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'};
for i=1:ncells_x            %For each sensory neuron,
    for j=1:ncells_y        %... and for each motor neuron,
        subplot(3,3,3*(i-1)+j)
        hold on             %... hold on to the graphics window,
        plot(z(:,i,j), [color{j} shape{i}]);  %... and plot the synaptic weight dynamics, with a chosen color and shape.
        legend(legendnames{3*(i-1)+j});
        hold off            %... then release the graphics window,
        [hfig(3*(i-1)+j)]=gca;
    end                     %... and close the for-loops.
end
xlabel('Time [steps]');  ylabel('Synaptic weights (z)')
%legend({'z_{11}'; 'z_{12}'; 'z_{13}'; 'z_{21}'; 'z_{22}'; 'z_{23}'; 'z_{31}'; 'z_{32}'; 'z_{33}'});
%linkaxes(hfig)
