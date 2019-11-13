function [xydata,violationcount]=cleanxy(xydata)
% basically makes sure timestamps arent repeated...
violations=find(diff(xydata(:,1))<=0);
xydata(violations,:)=[];
violationcount=numel(violations);