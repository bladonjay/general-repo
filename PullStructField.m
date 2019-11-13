function [fielddata] = PullStructField(struct,field)
%Get a field from a struct
% just because i wont remember
% JHB

fielddata=[struct(:).(field)];
end

