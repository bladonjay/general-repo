function spkmat = fixSPKMAT(samples, start,duration,varargin) 
% this function changes the sample window for a spkmat
% INPUTS
%   samples; sample matrix from spksevs file
%   start; the start time you want to use
%   duration; the duration of your sample
%   varargin; unit data
%OUTPUTS
%   spkmat: your new and improved spkmat

% originated from Sam Mckenzie, edited by John Bladon


nunits = 0;
spkmat = nan(size(samples,1),0);
for d = 1:nargin-1
    daysamples = samples(samples(:,10)==d,:);
    
    dayunitdata = varargin{d};
    if isfield(dayunitdata,'units')
        ndayunits = numel(dayunitdata.units);
    else
        ndayunits = numel(dayunitdata.neurons);
    end
    
    
    dayspkmat = nan(size(daysamples,1),ndayunits);
    for n = 1:ndayunits;
        
        
        
        if isfield(dayunitdata,'units')
            spk = dayunitdata.units(n).ts(:);
        else
            spk = dayunitdata.neurons{n}.timestamps(:);
        end
        
        
        for s = 1:size(daysamples,1)
         
                thissample = (spk(:,1)>=(daysamples(s,1) +start) & spk(:,1)<daysamples(s,1)+duration)/(duration);
            
            nspks = sum(thissample); % Count the Number of SPiKeS
            dayspkmat(s,n) = nspks; % Store the rate of spikes in the spike matrix.
        end
    end
    if(nunits < ndayunits)
        nunits = ndayunits;
        spkmat(:,end+1:nunits) = 0;
    end
    spkmat(samples(:,10)==d,1:ndayunits) = dayspkmat;

end