function [mergedstruct] = MergeStructs(newstruct,oldstruct,mergemethod)
% function [mergedstruct] = MergeStructs(newstruct,oldstruct,mergemethod)
% this function merges two structs.  Basically it will either append the
% new struct onto the bottom of the old, or it will update all the fields
% in the old struct to become those in the new

switch mergemethod
    case 'merge'
        % in this case, we will add to each field
        temp=[];
        
    case 'overwrite'
        % in this case we will overwrite any fields that are there, or add
        % new fields in, but keep old fields
        temp=[];
        
        
    case 'append'
        % in this case we will make sure the two structs have the same
        % fields, then we will tack the new struct onto the old struct
        temp=[];
        
        
end






end



function newstruct = checkfieldtypes(newstruct, refstruct)
fields = fieldnames(refstruct);

for j = 1:size(fields,1)
    % if new struct doesnt have the field, or its class is
    % wrong, new struct takes old structs field
    if(~isfield(newstruct,fields{j}) || ...
            ~isa(newstruct.(fields{j}),class(refstruct.(fields{j}))))
        newstruct.(fields{j}) = refstruct.(fields{j});
        % if newstruct has a field but its empty, or it has a nan
    elseif(isnumeric(newstruct.(fields{j})) && ...
            (isempty(newstruct.(fields{j})) || ...
            any(isnan(newstruct.(fields{j})))))
        newstruct.(fields{j}) = refstruct.(fields{j});
        warning(['Invalid number entered for ' fields{j} ', resetting to default value.'],'Invalid Entry');
    elseif(isstruct(newstruct.(fields{j})))
        newstruct.(fields{j}) = checkfieldtypes(newstruct.(fields{j}), refstruct.(fields{j}));
    end
end
end

