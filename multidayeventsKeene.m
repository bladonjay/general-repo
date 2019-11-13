function samples = multidayeventsKeene(varargin)
% Call this function with a list of cinedata variables, one for each day.
% samples = multidayevents(cinedata1, cinedata2, cinedata3, ...)
% to specify 'before' and 'after' times, add numers into the input
% arguments... for example:
% samples = multidayevents(before, cinedata1, cinedata2, ...) must be
% inputed in the order in which they were recorded
% samples = multidayevents(before, after, cinedata1, cinedata2, ...)
% samples = multidayevents([before, after], cinedata1, cinedata2, ...)

args = varargin; % Copy the list of input arguments.
tokeep = false(size(args)); % Keep a record of the inputs that are not structures.
for ii = 1:numel(args) % For each argument...
     if(isstruct(args{ii}))
      tokeep(ii) = true; % Record that this one is not a structure, and should be removed.
     end
    if(isnumeric(args{ii})) % Check if it is numerical (therefore not a structure)
      
        if(isscalar(args{ii})) % Check if it is one value.
            % If this is the first value we have gotten, call it 'before'.
            if(exist('before','var')~=1); before = args{ii};
                
            % If this is the second value we have gotten, call it 'after'.
            elseif(exist('after','var')~=1); after = args{ii};
                
            % If this is the third (or more), then give an error.
            else error('Too many numerical inputs.');
            end
        elseif(numel(args{ii})==2); % If provided as a vector
            % First value is 'before', second value is 'after'.
            if(exist('before','var')~=1 && exist('after','var')~=1);
                before = args{ii}(1);
                after = args{ii}(2);
            end
        end
    end
end
args = args(tokeep);

% If 'before' and/or 'after' were not included above, default to empty.
if(exist('before','var')~=1); before = []; end
if(exist('after','var')~=1); after = []; end


samples = zeros(0,15);

if any(tokeep)
for ii = 1:sum(tokeep);
    % andreaeventsKeene is poorly annotated and confusing, aviod usage
    samples = [samples; andreaeventsKeene(args{ii}, ii, before, after)]; %#ok<AGROW>
end
end
end