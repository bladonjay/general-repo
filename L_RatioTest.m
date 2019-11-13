function [IsolationStats] = L_RatioTest(mywfs,otherwfs,df)
% L_ratio = L_RatioTest(mywfs,otherwfs,varargin)
% runs an L-ratio test on my wfs vs other wfs using a bunch of candidate
% parameters

% Usage:
% INPUTS
%   mywfs is a m wfs by n components per wave for four tetrodes
%       if you have a dead channel, nan it out so you dont use the data but keep
%       the place holders in
%   otherwfs: same exact thing
%   df: trodality of this tetrode (default is 4
%   
% OUTPUTS:
% Struct containing:
%   IsoDist
%   L_Ratio
%   CovComp
%   CovAll

%  JHB

% not gonna use the parser for now
if (~exist('df','var') || size(temp)~=1)
    df=4;
end

%useparams=p.Results.params;

% i'm just gonna do it on all the params

% maybe make it a m waves by n pp wave by 4 waves 

allwfs=[mywfs; otherwfs];
wfdes=[ones(size(mywfs,1),1);zeros(size(otherwfs,1),1)];
for i=1:size(allwfs,2)
    if isnan(nanvar(mywfs(:,i)))
        kill(i)=1;
    else
        kill(i)=0;
    end
end
allwfs(:,kill==1)=[];

% now get the distance

% first get mahalanobis distance within
wmdist=sort(mahal(allwfs(wfdes==1,:),allwfs(wfdes==1,:)));

% now get mahalanobis distance noise
nmdist=sort(mahal(allwfs(wfdes==0,:),allwfs(wfdes==1,:)));
   
% isolation distance is the nth closest noise spike, where n is the number
% of spikes in that wave
IsoDist=nmdist(min([length(wmdist) length(nmdist)]));

% L ratio is the area inside the curve of i really dont know but below is
% the math
dfreedom=df; % should figure this out but its usually trodality
CDFnoise=cdf('Chisquare',nmdist,dfreedom);
if CDFnoise==0
    L_Ratio=nan;
else
    L_Ratio=nansum(1-CDFnoise)/length(nmdist);
end
%and now a covariance for the components
CovComp=cov(allwfs(wfdes==1,:));
% and cov for everything
CovAll=cov(allwfs);

IsolationStats.IsoDist=IsoDist;
IsolationStats.L_Ratio=L_Ratio;
IsolationStats.CovComp=CovComp;
IsolationStats.CovAll=CovAll;


end


