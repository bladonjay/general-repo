function [results,legend]=EventAnova(units,samples,cols,varargin)
% function [results]=EventAnova(units,samples,cols,varargin)
% A very general anova, where you pick your spikerates around the event,
% normalize them in some way, and then run an anova on it


p=inputParser;
% sampling window
windowvalid= @(a) numel(a)==2 & diff(a)>0;
defaultwin=[-1 2];
addOptional(p,'eventwindow',defaultwin); %windowvalid)

% transform
transformvalid= @(a) ischar(a);
addOptional(p,'transform','bymax');
% lets say the two options are zscore or normalize by max

parse(p,varargin{:});

window=p.Results.eventwindow;
transform=p.Results.transform;

% build the spike matrix
ratemat=eventspikematrix(samples(:,1),units,window(1),window(2));


if strcmpi(transform,'bymax')
    spkmat=NormalizeByMax(ratemat);
elseif strcmpi(transform,'zscore')
    spkmat=zscore(ratemat);
else
    spkmat=ratemat;
end

switch length(cols)
    case 0
        legend={'nothing to compute'};
        results=[];
        fprintf('Nothing to compute \n');
        return
    case 1
        legend={'A'};
    case 2
        legend={'A';'B';'AxB'};
    case 3
        legend={'A';'B';'C';'AxB';'AxC';'BxC';'AxBxC'};
    case 4
        legend={'A';'B';'C';'D','AxB';'AxC';'AxD';'BxC';'BxD';...
            'CxD';'AxBxC';'AxBxD';'AxCxD';'BxCxD';'AxBxCxD'};
end 

% now run stats on each column (unit)
for i=1:size(spkmat,2)
    % run anova on all the event types
    try
        [a,b,c]=anovan(spkmat(:,i),samples(:,cols),'model','full','sstype',2);
        %
        kill
        results(:,i)=a;
    catch
        fprintf('anova for cell %s doesnt work \n',strtrim(units.units(i).units));
        results(:,i)=nan;
    end
end
   


end

% First column is timestamp

% Second column is correct (1) / incorrect (0) (pot)
% Third column is left (1) / right (2)
% Fourth column is item set (context pair) (A/B = 1, C/D = 2)
% Fifth column is west (context 1) or east (context 2)
% Sixth column is position (1-4)
% Seventh column is Odor Pair (A/C = 1, B/D = 2)
% Eighth column is Odor (A = 1, B = 2, C = 3, D = 4)
% Ninth column is duration of the sample.
% Tenth column is the day number
% Eleventh column is actual trial #
% Twelth column is the # of samples on that pot
% Thirteenth column is rat being correct (dug on cor. no dig on incorrect)
% Fourteenth column is last sample of that pot for that trial)
% Fifteenth column is whether to exclude the trial


