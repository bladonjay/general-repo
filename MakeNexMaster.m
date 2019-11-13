function [nex,source] = MakeNexMaster(dirName)
% Builds a master NEX file for S Mackenzies nex2CMB script
% function[nex] = MakeNexMaster(dirName)
% input can be blank, or a directory name, output is a struct that can be
% read by NEX2CMBnew
% calls functions getALLNEXFiles and readNexFileM
% JHB 6-10-14


% if you didnt choose input, allow you to pick folder
if ~exist('dirName','var')
    dirName=uigetdir;
end


nexstruct=[];
fileList=getAllNEXFiles(dirName);
source=dirName;


%get all the nex data structs first
for i=1:length(fileList)
    disp(['Calculating' fileList{i}]);
    
    %read the nex files and put them into their individual fields
    [nex] = readNexFileM(fileList{i});
    nexstruct(i).nex=nex;
end

% Now combine them into a single nex struct
% cycle through all individual structs and add to existing fields in output
% struct
for i=1:length(nexstruct)
    %create the struct first
    if i==1
        nex=nexstruct(i).nex;
    
    % now add to the struct only if there are fields to add    
    else
        %make sure the end time is correct
        if nex.tend <  nexstruct(i).nex.tend
            nex.tend=nexstruct(i).nex.tend;
        end
        
        %if this is the behavioral data
        if isfield(nexstruct(i).nex,'events')
        
            % if the struct doesnt already have the data
            if ~(isfield(nex,'events'))

                % add the fields
                nex.events=nexstruct(i).nex.events;
                nex.markers=nexstruct(i).nex.markers;
            end
        end
            
        % if this file has neurons
        if isfield(nexstruct(i).nex,'neurons');
            % if its the first neuron file, create the field
            if ~(isfield(nex,'neurons'))
                nex.neurons=nexstruct(i).nex.neurons;
                
            % if its not the first neuron file, add to the field
            else
                nex.neurons=[nex.neurons;nexstruct(i).nex.neurons];
            end
        end
    end
end
        
        
end

