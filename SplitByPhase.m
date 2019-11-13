function [newunitdata] = SplitByPhase(unitdata,LFPts,LFPphases,bins)
%splits the unitdata struct up into more units based on theta phase
it=1;
if ~exist('bins','var')
    bins=6;
end
% bin the theta into bits
thetaparts=linspace(-pi,pi,bins+1);
% for each unit, grab its ts, get their phases
for i=1:length(unitdata)
    spikephases=interp1(LFPts,LFPphases,unitdata(i).ts);
    % and split into new units by phase
    for j=1:max(bins)
        newunitdata(it)=unitdata(i);
        newunitdata(it).ts=unitdata(i).ts(spikephases>thetaparts(j) & spikephases<=thetaparts(j+1));
        it=it+1;
    end
end

end
