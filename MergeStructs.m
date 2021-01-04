function [mergedstruct] = MergeStructs(oldstruct,newstruct,mergemethod)
% function [mergedstruct] = MergeStructs(newstruct,oldstruct,mergemethod)
% this function merges two structs.  Basically it will either append the
% new struct onto the bottom of the old, or it will update all the fields
% in the old struct to become those in the new. The resulting struct will
% have the fields of the new struct.  To merge two structs completely, you
% can call this function twice, reversing the order fo the structs.
% old strcut can be an n element struct, but new has to be 1


% update makes a new struct where all fields are present, but overlapping
% fields take the data from the newstruct

% addoverlap just tacks onto the old struct a new row for the new struct,
% and it only takes the fields shared by both structs

% addOne just takes the old fields, adds a new row from the new struct, and
% leaves new fields not part of the old struct alone.  if the new struct
% doesnt hvae old fields, they remain blank

% addAll takes all fields from all structs and tacks on.  here you will get
% empty fields in the old and new structs if they dont ahve those fields.

if ~exist('mergemethod','var')
    mergemethod='addon';
end

mergemethod=lower(mergemethod);
oldlength=length(oldstruct);

newfields=fieldnames(newstruct);
oldfields=fieldnames(oldstruct);
allfields=unique([newfields; oldfields]);
sharedfields=intersect(newfields,oldfields);
oldonly=setdiff(oldfields,newfields);
newonly=setdiff(newfields,oldfields);
switch mergemethod
    % in this case, all fields that exist in the first struct are put in the
    % second struct, if they dont exist in the first struct, the second are
    % kept
    case 'update'
        for i=1:length(allfields)
           if exist(oldstruct,(allfields{i}))
               if ~isempty(oldstruct.(allfields{i}))
                   mergedstruct.(allfields{i})=oldstruct.allfields{i};
               else
                   mergedstruct.(allfields{i})=newstruct.allfields{i};
               end
           else
               mergedstruct.(allfields{i})=newstruct.allfields{i};
           end
        end
        % in this case, we build a 2 row struct, hwere only the shared
        % features are kept
    case 'addoverlap'
        for i=1:length(sharedfields)
            for j=1:length(oldstruct)
                mergedstruct(j).(sharedfields{i})=oldstruct(j).(sharedfields{i});
            end
            for j=1:length(newstruct)
                mergedstruct(end+1).(sharedfields{i})=newstruct.(sharedfields{i});  
            end
        end
        % in this case this just adds a new struct on bottom but only uses
        % the old fields
    case 'addon'
        for i=1:length(oldfields)
            for j=1:length(oldstruct)
                mergedstruct(j).(oldfields{i})=oldstruct(j).(oldfields{i});
            end
            for j=1:length(newstruct)
                if isfield(newstruct,oldfields{i})
                    mergedstruct(oldlength+j).(oldfields{i})=newstruct(j).(oldfields{i});
                else
                    mergedstruct(oldlength+j).(oldfields{i})=[];
                end
            end
        end
        % in this case we just add alln of new struct to old struct
    case 'addall'
        % first add all
        for i=1:length(oldfields)
            for j=1:length(oldstruct)
                mergedstruct(j).(oldfields{i})=oldstruct.(oldfields{i});
            end
            for j=1:length(newstruct)
                if isfield(newstruct,oldfields(i))
                    mergedstruct(oldlength+j).(oldfields{i})=newstruct(j).(oldfields{i});
                else
                    mergedstruct(oldlength+j).(oldfields{i})=[];
                end
            end
        end

        for i=1:length(newonly)
            for j=1:length(oldstruct)
                mergedstruct(j).(newonly{i})=[];
            end
            for j=1:length(newstruct)
                mergedstruct(oldlength+j).(newonly{i})=newstruct(j).(newonly{i});
            end
        end
        
end

mergedstruct=orderfields(mergedstruct);





end




