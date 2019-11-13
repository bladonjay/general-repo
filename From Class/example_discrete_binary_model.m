function Sj = example_discrete_binary_model(S1, S2)

  % Set the weights.
  w1j= 0.5;
  w2j=-0.5;
  
  % Define the activity at neuron j, xj,
  xj = S1*w1j + S2*w2j;
  
  % Apply the binary threshold,
  if xj > 0
      Sj = 1;
  else
      Sj = 0;
  end
  
end