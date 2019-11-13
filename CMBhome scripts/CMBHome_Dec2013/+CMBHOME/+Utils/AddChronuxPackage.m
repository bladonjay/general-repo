function AddChronuxPackage

if ~isempty(findstr('chronux', path)), return; end % check that it isnt already there

a = what('CMBHOME');

addpath(path, genpath(fullfile(a.path(1:end-8), 'chronux'))); % adds chronux package