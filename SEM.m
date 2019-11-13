function [ster,n] = SEM (X,dim)
if ~exist('dim','var')
    [~,dim]=max(size(X));
end
n=size(X,dim);
n_nan=sum(isnan(X),dim);
n=n-n_nan;
dv=nanstd(X,[],dim);
ster=dv./sqrt(n);
end