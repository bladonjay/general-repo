%%  NE204 - Lab 5:  The integrate and fire neuron.
%   In this lab, we will study the integrate and fire (I&F) model neuron.  We
%   will simulate this model in MATLAB, and get a sense for how it behaves.
%   We'll also manipulate the model a bit, and gain some intuition for how
%   the model changes, and the biological reasons for these changes.

%%  Preliminaries.
%   Text preceded by a '%' indicates a 'comment'.  This text should appear
%   green on the screen.  I will use comments to explain what we're doing 
%   and to ask you questions.  Also, comments are useful in your own code
%   to note what you've done (so it makes sense when you return to the code
%   in the future).  It's a good habit to *always* comment your code.

%%  Part 1:  Preliminaries - Get the I&F code  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Download from Blackboard the following code,
%
%     leaky_integrate_fire.m
%
%   Save this code in a place you can locate, and then navigate MATLAB to
%   this place.
%
%   Once you're there, execute the following code:

leaky_integrate_fire

%    Doing so, a beautiful graphical user interface (GUI) should appear.
%    We'll use this interface to study the I&F model in this lab.

%IN LAB Q:  Click around on this GUI and explore.  What can you make it do?

%%  Part 2:  I&F experiments Round 1: Vary the constant current injected.
%   Now that we've lauched the GUI, let's perform some experiments.  To
%   start, let's increase the constant current injected to the model.  

%IN LAB Q:  Begin with the default settings.  Here, I=0.  How do you
%*expect* the voltage dynamics to evolve in time?  Run the simulation, and 
%examine how the voltage dynamics (V) evolve in time.  What do you observe?
%Does it match your expectations?

%IN LAB Q:  What is the *rest* voltage in this model?

%IN LAB Q:  Increase the constant current input I, and examine how the
%dynamics of V (the voltage) change.  What do you observe?  Can you make
%the model "spike"?  What do the spikes "look like" in the LIF model?

%IN LAB Q:  Determine how much constant input current I must be injected to
%produce an action potential (or "spike") in the model.

%IN LAB Q:  What is the *threshold voltage* in this model?

%IN LAB Q:  What is the *reset voltage* in this model?

%IN LAB Q:  Examine the current flow through the resistor (I_res in red) and
%through the capacitor (I_cap in blue).  Describe how these two quantities
%change leading up to a spike. 

%%  Part 3:  I&F experiments Round 2:  Vary current pulses and widths.
%   In Part 2, we varied the constant current injected.  Let's now consider
%   how the I&F model responds to pulses of injected current.  To do so,
%   select the button,
%
%   "current pulse:"
%
%   In the "input current I(t)" window.

%IN LAB Q:  What do T and I represent? (i.e., What does pulse width and
%pulse height mean?)

%IN LAB Q:  Simulate the model neuron with current pulses of different widths
%(T) and heights (I) such that the product T*I is constant (18 ms*mA). Is it
%easier to get the neuron to fire with (a) a brief, intense pulse or (b) an
%extended, low amplitude pulse? (Hint: check these cases by switching two
%integers that multiply to 18 ms*mA, keeping in mind that the program only
%allows input currents up to 10 mA.)

%IN LAB Q:  Explain, intuitively, why changing the shape of a current pulse
%matters for a leaky integrate and fire neuron in terms of cell
%properties. (HINT:  Think about the target voltage.)

%% Part 4:  A very simple I&F model.
%  In class, we wrote down the differential equation for the LIF model, as
%  well as the threshold & reset condition.  So far, we've used the fancy
%  GUI software to compute the solution.  That's nice because it makes a great
%  visualization, and we can get a sense for how the system behaves by
%  quickly adjusting parameters and seeing what happens.
%
%  But, how is the solution actually computed?  This model is not too
%  difficult to simulate "by hand":  we can write our own code to simulate
%  the LIF model which, as we'll see, involves a for-loop and an
%  if-statement.  From Blackboard download the function,
%
%  LIF.m
%
%  In lab, we'll discuss this function, and get a sense for how it behaves.
%  The same simulation strategy will apply to other more complicated models
%  we will consider in this course.  So, it's useful to get some intuition
%  now, for this relatively simple model.
%
%  Let's simulate the LIF model and see how it behaves.  To start,
%  consider,

Iinput=3;
R=10;
C=0.6;
Vrest=-70;
Vthreshold=-50;
Vreset=-80;

[V,t,nspikes] = LIF(Iinput,R,C,Vrest,Vthreshold,Vreset);

%  Notice, that the function LIF takes 6 inputs.  We set the values for these
%  six inputs in the lines above the function call.  Also notice that the
%  function LIF returns 3 output:
%
%  V = the LIF voltage.
%  t = the time axis (useful for plotting).
%  nspikes = the number of spikes observed.
%
%  To display the results, let's plot the voltage versus time,

plot(t,V)

%  For these parameter settings, the LIF neuron generates lots of spikes.
%  To see how many, let's look at the third output of the LIF function,

nspikes

%IN LAB Q:  Compare the value of 'nspikes' to the plot of V versus t.  Do
%the number of spikes match?  (HINT:  They should!)


%% Part 4: Exploring our toolbox of equations.

% We have a number of equations to allude to when we are trying to
% understand how our differential LIF equation works.

% lets start with Ohms Law: V=IR
%
% IN LAB Q:  Intuitively, what is V, what is I and what is R?
%
% IN LAB Q: Describe in the biological system what it would mean if we
% reduce R.


% Now lets look at our time constant: t=RC
%
% IN LAB Q: intuitively, lets describe each of these terms
%
% IN LAB Q: An easy analogy for this equation is to think of a river that splits into a
% narrow part and a resevoir with a dam.  Lets see if we can describe what happens when we
% INCREASE C in this equation


% Now lets look at *target voltage*   Vtarget=Iin*R + Vrest
%
% IN LAB Q: again, lets describe each of these variables
%
% Now lets think of what happens when we INCREASE R

%% Part 5: Dissecting our LIF function.

% lets enter the function LIF to see how it works
edit LIF


% Can you find the time constant and the target voltage in this
% function?




