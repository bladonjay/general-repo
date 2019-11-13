%Implement the Hodgkin-Huxley equations.
%See Gerstner and Kistler, Spiking Neuron Models, 2002, Section 2.2.
%You'll see I've scaled the voltage by 65 in the equation that updates V
%and the auxillary functions.  Hodgkin and Huxley set the resting voltage
%of the neuron to 0 mV and we've set it here to -65 mV (the value accepted
%today).

%INPUTS
%    I0 = input current.
%    T0 = total time to simulate (in [ms]).
%
%OUTPUTS
%    V = the voltage of neuron.
%    m = activation variable for Na-current.
%    h = inactivation variable for Na-current.
%    n = activation variable for K-current.
%    t = the time axis of the simulation (useful for plotting).

function [V,m,h,n,t] = HH0(I0,T0)

  dt = 0.01;
  T  = ceil(T0/dt);
  gNa0 = 120;  %mS/cm^2 (sodium channel conductance)
  ENa  = 115;  %mV (equilibrium potential for sodium)
  gK0  = 36;   %mS/cm^2
  EK   = -12;  %mV
  gL0  = 0.3;  %mS/cm^2 % leak current
  EL   = 10.6; %mV % equilbrium pot for leak current

  % set up outputs
  t = (1:T)*dt;
  V = zeros(T,1);
  m = zeros(T,1);
  h = zeros(T,1);
  n = zeros(T,1);
  
  % initial conditions
  V(1)=-70.0; % start at -70 mV
  m(1)=0.05; % Sodium activation
  h(1)=0.54; % Sodum inactivation
  n(1)=0.34; % activation for K current (closes based on time)
  
  % now run our simulation
  for i=1:T-1 
      % voltage for next step
      V(i+1) = V(i) + dt*(gNa0*m(i)^3*h(i)*(ENa-(V(i)+65))...
                        + gK0*n(i)^4*(EK-(V(i)+65)) + ...
                          gL0*(EL-(V(i)+65)) + I0);
      
      % sodium channels open
      m(i+1) = m(i) + dt*(alphaM(V(i))*(1-m(i)) - betaM(V(i))*m(i));
      % how many sodium channels are inactivated?
      h(i+1) = h(i) + dt*(alphaH(V(i))*(1-h(i)) - betaH(V(i))*h(i));
      % how many potassium channels are open?
      n(i+1) = n(i) + dt*(alphaN(V(i))*(1-n(i)) - betaN(V(i))*n(i));
  end
  
end

%Below, define the AUXILIARY FUNCTIONS alpha & beta for each gating variable.

% the percent of sodium channels that are activated
function aM = alphaM(V)
aM = (2.5-0.1*(V+65)) ./ (exp(2.5-0.1*(V+65)) -1);
end

% percent of sodium channels NOT activated
function bM = betaM(V)
bM = 4*exp(-(V+65)/18);
end

% sodium inactivation state (cant be activated)
function aH = alphaH(V)
aH = 0.07*exp(-(V+65)/20);
end

% sodium not inactivation state
function bH = betaH(V)
bH = 1./(exp(3.0-0.1*(V+65))+1);
end

function aN = alphaN(V)
aN = (0.1-0.01*(V+65)) ./ (exp(1-0.1*(V+65)) -1);
end

function bN = betaN(V)
bN = 0.125*exp(-(V+65)/80);
end
