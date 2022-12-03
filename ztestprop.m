%ztestprop
%z test for two proportions
%x1 and x2 are the number of 'positive' results in each sample
%n1 and n2 are the sizes of the two samples
function [p,z] = ztestprop(x1, x2, n1, n2)

p1 = x1/n1;
p2 = x2/n2;

prob = (x1+x2)/(n1+n2);
pinv = 1-prob;

z = (p1-p2)/sqrt(prob*pinv*((1/n1)+(1/n2)));

p = normcdf(z);
if p == 0 || p == 1
    p = normcdf(z,'upper')
end
    


