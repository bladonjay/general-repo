function [pval,chi2stat,df]=chi2indep(datamat,df,expected)
% function [pvalue,chi2stat=chi2indep(datamat,df,expected)
%performs a chi square test of independence, basically whether the observed
% proportions are significantly different from chance.  E.G. are the
% columns and rows
% Inputs:
%   Datamat: The data matrix, data matrix of at least 2 columns and at least 2 rows
%   df: degrees of freedom, usually rows-1 * cols-1
%   Expected:  if you have expected proportions, the expected mat has to be the same
%   size as the datamat
% Outputs:
%   pval: the p value of the test
%   Chi2stat: the raw chi squared statistic thats compared against the
%   distribution


[rows,cols]=size(datamat);
if ~exist('df','var')
    df=(rows-1)*(cols-1);
end
% if no expected mat, calculate the expected
if ~exist('expected','var')
    % tabulate sums
    colsums=sum(datamat);
    rowsums=sum(datamat,2);
    N=sum(colsums);
    
    % now fill in each variable as if all proportions were even
    for i=1:rows
        for j=1:cols
            expected(i,j)=(rowsums(i)*colsums(j))/N;
        end
    end    
    
end
% now calculate sum((obs-exp)^2/exp)
chi2stat=sum(sum((datamat-expected).^2./expected));
% and compare to a chi square distribution
pval= 1.0000- chi2cdf(chi2stat,df);
end
