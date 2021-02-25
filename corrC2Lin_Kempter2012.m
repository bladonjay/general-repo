function [a,phi_knot,rho,pval,ev] = corrC2Lin_Kempter2012(x,phi)
% [a,phi_knot,rho_c,p] = utils.corrC2Lin_Kempter2012(x,phi)
%
% inputs: 
%   x - Linear variable
%   y - Angular variable (in radians, -pi:pi or 0:2pi)
% 
% outputs:
%   a        - Slope (rads / x unit)
%   phi_knot - Angular intercept
%   rho    - Correlation coefficient
%   p        - probability that slope is significantly non-zero
%   ev      -expected phase for each spike
% See R. Kempter 2012, Journal of Neuroscience Methods 
% wchapman 05.12.2012

%% Clean the data
bads = isnan(x) | isnan(phi);
x(bads) = [];
phi(bads) = [];

n = length(phi);

%% Fit the line by MSE
% seek the best slope by minimizing the slope and the offset
% this will get you the slope and offset

%params = fminsearch(@(params) costFunct(params, x, phi), [pi, pi / (max(x)-min(x))]); %used for original for 8s sessions analysis

% use this for diff tr speed precession slopes
% error: Exiting: Maximum number of function evaluations has been exceeded,
% current value 56.39
options=optimset('MaxFunEvals',10000,'MaxIter',10000); %added MaxIter for time vs distance phase rasters
params = fminsearch(@(params) costFunct(params, x, phi), [pi, pi / (max(x)-min(x))],options);

a = params(2); 
phi_knot = mod(params(1),2*pi);
ev=mod(x*a+phi_knot,2*pi);


%% Calculate correlation coefficient
% this is independent of slope and offset because its sin over cosine.
n = length(x);
rxs = corr(x,sin(phi));
rxc = corr(x,cos(phi));
rcs = corr(sin(phi),cos(phi));
% compute angular-linear correlation (equ. 27.47)
rho = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2));

% compute pvalue
pval = 1 - chi2cdf(n*rho^2,2); %


end

function cost = costFunct(params, x, phi)

    phi_hat = params(1) + x*params(2);
    
    c1 = (cos(phi_hat) - cos(phi)).^2;
    c2 = (sin(phi_hat) - sin(phi)).^2;
    
    cost = sum(sqrt(c1+c2));

end

function lv = lambda(i,j,x,phi,phi_bar,theta,theta_bar)

lv = (1/length(phi))*sum((sin(phi-phi_bar)) .^i .* (sin(theta - theta_bar)) .^j);

end