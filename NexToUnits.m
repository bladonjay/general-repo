function [units,totcel] = NexToUnits(dirName,includeWaves)
% function [Units] = NEXtoUnits(filePath,includeWaves)
% this function takes the datafiles generated on offlineSorter and
% combines them into a single "units" struct.  

if ~exist('includeWaves','var')
    includeWaves=0;
end

includeWaves=any(includeWaves);

nexFiles=dir([dirName '/**/*.nex']);
units=struct;
totcel=1;
for i=1:length(nexFiles)
    
    fprintf('Reading %s \n', nexFiles(i).name);
    datafile= readNexFileM(fullfile(nexFiles(i).folder,nexFiles(i).name));
    mytetstr=regexp(nexFiles(i).name,'.nt[0-9]+','match');
    mytet=str2double(mytetstr{1}(4:end));
    
    if isfield(datafile,'neurons')
        fprintf('found %d units \n', length(datafile.neurons));
        for j=1:length( datafile.neurons)
            
            % trim blank space before and after
            tempname= strtrim(datafile.neurons{j}.name);
            tempname(tempname=='_')='-'; % underscore is an ascii special character
            % remove all the spaces
            units(totcel).name = tempname(tempname~='_');
            units(totcel).tet = mytet;
            units(totcel).cellnum=j;
            % get timestamps
            units(totcel).ts = datafile.neurons{j}.timestamps;
            % pull in a waveform too
            if isfield(datafile,'waves') && includeWaves==1
                units(totcel).waves=datafile.waves{j}.waveforms;
            end
            totcel=totcel+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% interneuron filter goes her %%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if you dont want to pull in all of the waveforms, you would
            % calculate waveform stats (isodist, L ratio, mean wave shape)
        end
    end


end

