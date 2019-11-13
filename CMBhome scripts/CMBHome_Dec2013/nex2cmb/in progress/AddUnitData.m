function [spikedata, nex] = AddUnitData(nexfile,timestamps)
%Creates a spike data object in the CMBHOME format
% JHB
% I think i need to add a variable dataobject in ( could probably get this
% from a file, or if built into a wrapper)

eventstokeep = {};
 batch = false;
 nlx = false;
t=timestamps;
% %refers to Spk_i?
% 
% for a = 1:numel(varargin)
%     if(iscellstr(varargin{a}))
%         % 'eventstokeep' dictates which events from the NEX file are copied into
%         % the CMBHOME object. If empty, all events are kept. To exclude all events,
%         % include one element in 'eventstokeep' that won't match any other event,
%         % for example: eventstokeep = {''};
%         eventstokeep = lower(varargin{a});
%     elseif(islogical(varargin{a}) && isscalar(varargin{a}))
%         % Batch mode, don't display anything
%         batch = varargin{a};
%     elseif(ischar(varargin{a}) && any(strcmpi(varargin{a},{'nvt','nlx'})))
%         % This is for a neuralynx file I believe
%         nlx = true;
%     else error('NEX2CMB:UnrecognizedInput', 'Argument %d unrecognized.',a+1);
%     end
% end

% If no input file name is provided then pass an empty string to next level
if(nargin < 1 || isempty(nexfile)); nexfile = ''; end

% Read the NEX file into a structure.
if(isstruct(nexfile));
    nex = nexfile;
    if(isfield(nex,'nexfileinfo') && numel(nex.nexfileinfo)>=1 && ...
            isfield(nex.nexfileinfo,'name'))
        if(numel(nex.nexfileinfo)>1)
            names = char({nex.nexfileinfo(:).name}');
            samechar = all(diff(names+0,1,1)==0,1);
            nexfile = names(1,samechar);
            if(isempty(nexfile)); nexfile = nex.nexfileinfo(1).name; end
        else
            nexfile = nex.nexfileinfo.name;
        end
        if(~batch); fprintf(1,'Found NEX file name: %s\n',nexfile); end
    else
        nexfile = inputname(1);
        if(~batch); fprintf(1,'Using NEX file name: %s\n',nexfile); end
    end
else nex = readNexFileM(nexfile);
end


if isfield(nex,'neurons')
    names = cellfun(@(x) {x.name}, nex.neurons);
    
    % Get the tetrode numbers for each neuron.
    % This assumes that the tetrode number is the first block of consecutive
    % numbers within the tetrode name.
    tetnum = str2double(regexp(names,'([0-9])+','once','match'));
    
    % Get the unit number for each neuron.
    % This assumes that the units are named such that the last letter of the
    % name denotes the unit number: a = 1, b = 2, c = 3, etc. This is how
    % Offline Sorter names units when using the 'Export to NEX' option.
    unitnum = cellfun(@(x) x(end)-96,names);
    
    
    %cineplex studio does not give tetrodes so delete repeats
    units = unique([tetnum unitnum],'rows');
    
    for jj = 1:size(units,1)
       keep(jj) = find(ismember( [tetnum unitnum],units(jj,:),'rows'),1,'first');
        
    end
    
    nex.neurons = nex.neurons(keep);
    tetnum = tetnum(keep);
    unitnum = unitnum(keep);
    for ii = 1:size(nex.neurons,1)
        if(~batch); fprintf(1,'%s is tetrode %d unit %d\n',names{ii},tetnum(ii),unitnum(ii)); end
        spk_ts = nex.neurons{ii}.timestamps;
        spk_i = time2ind(t,spk_ts);
        spk_ts = spk_ts(spk_i>0);
        spk_i = spk_i(spk_i>0);
        %placeholder for spk_i
        %spk_i=1;
        
        spikedata(ii)= CMBHOME.Spike('i',spk_i,'ts',spk_ts,'tet_name',names{ii});
    end
end



end

