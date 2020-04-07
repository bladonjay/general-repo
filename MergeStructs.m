function [mergedstruct] = MergeStructs(oldstruct,newstruct,mergemethod)
% function [mergedstruct] = MergeStructs(newstruct,oldstruct,mergemethod)
% this function merges two structs.  Basically it will either append the
% new struct onto the bottom of the old, or it will update all the fields
% in the old struct to become those in the new. The resulting struct will
% have the fields of the new struct.  To merge two structs completely, you
% can call this function twice, reversing the order fo the structs.

newfields=fieldnames(newstruct);
oldfields=fieldnames(oldstruct);
allfields=unique([newfields; oldfields]);
switch mergemethod
    % in this case, all fields that exist in the old struct are puit in the
    % new struct
    case 'keep'
        for i=1:length(allfields)
           if exist(oldstruct,(allfields(i)))
               if ~isempty(oldstruct.(allfields(i)))
                   mergedstruct.(allfields(i))=oldstruct.allfields(i);
               else
                   mergedstruct.(allfields(i))=newstruct.allfields(i);
               end
           else
               mergedstruct.(allfields(i))=newstruct.allfields(i);
           end
        end
        % in this case, the new struct overwrites the old struct, but if
        % the new struct doesnt have those fields, they're imported from
        % the old struct
    case 'update'
        
        temp=[];

end






end




