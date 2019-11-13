function [Session] = NexToSession(filepath,inputsession)
% function [Session] = NexToSession(datafiles)
% this function takes the datafiles generated on Cineplex and Sorter and
% combines them into a single "session" struct.  It will not overwrite data, but it
% will combine all into a single session.  It will take the name of the
% first next file.
% its probably not a good idea to cat the sessions that come out of this
% because they will not have congruent fields

% the fields you'll get in the session are as follows:
% coords
% neurons
% rawevents
if exist('inputsession','var')
    Session=inputsession;
else
    Session=struct('coords',[],'unitdata',[],'rawevents',[]);
end

for i=1:length(filepath)
    % grab the data
    datafile{i} = readNexFileM(filepath{i});
    fprintf('Reading %s \n \n', filepath{i}(find(filepath{i}=='\',1,'last')+1:end-4));
   
    % there will be three types of data in these files that we care about,
    % spike data, coordinate data, and event data.
    
    % 1. Coordinate data, do you want it and is it already there?
    if  isfield(datafile{i},'markers')
        
            % if this file has markers, it may be the tracking coords
            markernames=cellfun(@(a) a.name,datafile{i}.markers,'uni',0);
            % Lets see what we can nab
            if any(cell2mat(strfind(markernames,'DVT')))
                fprintf('Found DVT coords \n');
                if ~isfield(Session.coords,'DVTcoords')
                    fprintf('Adding DVT coords \n');
                    Session.coords.DVTcoords=readNexCoords(datafile{i});
                else, Session.coords(end+1).DVTcoords=readNexCoords(datafile{i});
                end
            elseif any(cell2mat(strfind(markernames,'AVI')))
                fprintf('found AVI coords \n');
                if ~isfield(Session.coords,'AVIcoords')
                    fprintf('Adding AVI coords \n');
                    Session.coords.AVIcoords=readNexCoords(datafile{i});
                else, Session.coords(end+1).AVIcoords=readNexCoords(datafile{i});
                end
            elseif any(cell2mat(strfind(markernames,'PLX')))
                fprintf('found PLX coords \n');
                if ~isfield(Session.coords,'PLXcoords')
                    fprintf('Adding PLX coords \n');
                    Session.coords.PLXcoords=readNexCoords(datafile{i});
                else,  Session.coords(end+1).PLXcoords=readNexCoords(datafile{i});
                end
            else
                fprintf('Tracking marker name ''%s'' not recognized \n', markernames{1});
            end
            % do strobe separate, its not mutually exclusive
            if any(cell2mat(strfind(markernames,'Strobe')))
                fprintf('found Strobe coords \n');
                % if we have any coords
                if isfield(Session,'coords') 
                    % if we have coords but no strobe
                    if ~isfield(Session.coords,'Strobe')
                        fprintf('Adding Strobe coords \n');
                        Session.coords.strobe=readNexCoords(datafile{i});
                    end
                else
                    Session.coords.strobe=readNexCoords(datafile{i});
                end
            else
                fprintf('Found no separate Strobe Marker \n');
            end
        
    end
    
    % 2. now look for spike data
    % if there are neurons, or units, there should be only one of these
    % files, but we can combine need be later
    if isfield(datafile{i},'neurons')
        fprintf('found %d units \n', length(datafile{i}.neurons));
        if isfield(Session.unitdata, 'units')
            fprintf('already have %d neurons \n',length(Session.unitdata.units));
            oldunits=length(Session.unitdata.units);
        else
            oldunits=0;
        end
        for j=1:length( datafile{i}.neurons)
            % trim blank space before and after
            tempname= strtrim(datafile{i}.neurons{j}.name);
            % remove all the spaces
            Session.unitdata.units(j+oldunits).name = tempname(tempname~='_') ;
            % get timestamps
            Session.unitdata.units(j+oldunits).ts = datafile{i}.neurons{j}.timestamps;
            % pull in a waveform too
            if isfield(datafile{i},'waves')
                Session.unitdata.units(j+oldunits).wave=datafile{i}.waves{j}.waveforms;
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% interneuron filter goes her %%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
       
    end
    
    
    % 3. now look for flag data e.g. events, collect all events you can
    % find
    if isfield(datafile{i},'events')
        fprintf('found %d events \n', length(datafile{i}.events));
        for j=1:length(datafile{i}.events)
            if ~exist('rawevents','var')
                rawevents=datafile{i}.events{j};
            else
                rawevents(length(rawevents)+1)=datafile{i}.events{j};
            end
            % tag onto any found events
        end
        Session.rawevents=[Session.rawevents rawevents];
    end
  
    
    % 4. find LFPs and choose favorites
    if isfield(datafile{i},'contvars')
        
        
        contnames=cellfun(@(a) a.name,datafile{i}.contvars,'UniformOutput',false);
        choice=questdlg('Do you want to add all, some, or none of the contvars',...
            'Menu for LFP Data','All','Half','Choose/none','All');
        switch choice
            case 'Choose/none'
        [checked]=checkBox(contnames);
        % if any chosen, add all the chosen into the LFP
        if any(checked)
            for j=1:length(checked)
                Session.LFP(j).data=datafile{i}.contvars{checked(j)}.data;
                Session.LFP(j).name=datafile{i}.contvars{checked(j)}.name;
            end
        end
            case 'All'
                for j=1:length(contnames)
                    Session.LFP(j).data=datafile{i}.contvars{j}.data;
                    Session.LFP(j).name=datafile{i}.contvars{j}.name;
                end
            case 'Half'
                for j=1:floor(length(contnames)/2)
                    Session.LFP(j).data=datafile{i}.contvars{j*2}.data;
                    Session.LFP(j).name=datafile{i}.contvars{j*2}.name;
                end
        end
    end
    
end

% everything can just write over, datafile is the only one that needs to be
% cleared from workspace
clear datafile rawevents


end

