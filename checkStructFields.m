function newstruct = checkStructFields(newstruct, refstruct)
%function newstruct = checkStructFields(newstruct, refstruct)
% this function takes an old struct and updates all fields with the new struct
% if the new struct is missing fields in the old struct, they remain
% unchanged
% if the new struct has fields the old struct doesnt, they are created and
% added


            fields = fieldnames(refstruct);
            
            for j = 1:size(fields,1)
                % if new struct doesnt have the field, or its class is
                % wrong, new struct takes old structs field
                if(~isfield(newstruct,fields{j}) || ...
                        ~isa(newstruct.(fields{j}),class(refstruct.(fields{j}))))
                    newstruct.(fields{j}) = refstruct.(fields{j});
                % if newstruct has a field but its empty, or it has a nan
                % you keep the old struct field
                elseif(isnumeric(newstruct.(fields{j})) && ...
                        (isempty(newstruct.(fields{j})) || ...
                        any(isnan(newstruct.(fields{j})))))
                    newstruct.(fields{j}) = refstruct.(fields{j});
                    warning(['Invalid number entered for ' fields{j} ', resetting to default value.'],'Invalid Entry');
                % finally, if the field in the struct is a struct,
                % recursively run this function
                elseif(isstruct(newstruct.(fields{j})))
                    newstruct.(fields{j}) = checkStructFields(newstruct.(fields{j}), refstruct.(fields{j}));
                end
            end
        end
end

