function [fileList] = getAllFiles(dirName,extension)

  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  %fileData = dirData(~dirIndex);
  
  if exist('extension','var')
      if isnumeric(extension)
          fileList=fileList(cellfun(@any,strfind(fileList,'.mat'))); %only take matlab files
      else
          try
              fileList=fileList(cellfun(@any,strfind(fileList,extension)));
          catch
              warning('couldnt use extension');
          end
      end
  end

  if ~isempty(fileList)
      fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
          fileList,'UniformOutput',false);
  end
  
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
  %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
      nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
      if exist('extension','var')
      fileList = [fileList; getAllFiles(nextDir,extension)];  %# Recursively call getAllFiles
      else
          fileList = [fileList; getAllFiles(nextDir)];  %# Recursively call getAllFiles
      end
  end
  
  end