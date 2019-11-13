function [dash] = WhichDash
%dash=WhichDash creates a dash variable so you dont have to code a few extra lines of code
% works for macs and pcs
if ismac
    dash = '/';
elseif ispc
    dash = '\';
end

end

