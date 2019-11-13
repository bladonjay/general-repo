function Y  = fisherTransform(X)
X(X==1) = 0.999;
X(X==-1) = -0.999;
Y = .5 * log((1+X)./(1-X));



end