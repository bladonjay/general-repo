function d = dprime(X,varargin)
if ~isempty(varargin)
    Y = varargin{1};
    
    d = abs(nanmean(X) - nanmean(Y))./sqrt(.5*(nanvar(X)+nanvar(Y)));
else
    d = (nanmean(X) )/sqrt(nanvar(X));
    
    %d = 1-nanmean(Y)/nanmean(X);
end