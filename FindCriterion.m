function [CritTrial,CritIndex] = FindCriterion(trialnumber,performance,varargin)
%U function [Critical trial] = FindCriterion(trialnumber,performance)
% finds which trial the rat reaches criterion in a given task
% INPUTS
%   trialnumber: a vector of numbers, 1:number of trials
%   performance: a vector of 1s and 0's the same length as trials,
%                1 is correct 0 is incorrect
%   Optionals: 'window', default 20, how many trials is your window, has to
%                   be even
%              'crit', default 70% what percentage is criterion, use
%                   decimal (.7)
% OUTPUTS
%   CritTrial: the trial number where the rat hit criterion
%   Critindex: the index in your vector of trials where the criterion was
%   hit


% criterion is the midpoint over the trial window that the rat got correct

p=inputParser;
addOptional(p,'window',20);
addOptional(p,'crit',.7);
parse(p,varargin{:});

criterion=p.Results.crit;
window=p.Results.window;
% find the last sample for that trial (final choice)
[~,trialindices]=unique(trialnumber,'last');

lastsample=performance(trialindices);

% from trial 11:end
for q=window/2+1 : length(lastsample)-window/2
    % get percentage of those trials
    percentage=sum(lastsample(q-window/2 : q+window/2))/window;
    if percentage>criterion
        CritTrial=q; CritIndex=trialindices(q);
        break
    end
end


% now check our criterion to see if its missed or bad
if ~exist('CritTrial','var')
    CritTrial=nan; CritIndex=nan;
end

end

