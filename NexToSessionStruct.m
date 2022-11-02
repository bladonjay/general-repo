function [Session, unitdata, flags, coords, datafile] =NexToSessionStruct(dirName,varargin)
%[Session, units, flags, coords] =NexToSessionStruct(dirName,varargin);
% makes Session struct from a folder of nex files, see below
% this is for the Delay Task

%% Annotations

% INPUTS:
%   dirName; path of directory full of one session of nex files.  This is
%   usually one flags file and one units file
%

% OUTPUTS:
%   Session: struct with Coords, fixed Coords, flags, samples, responses,
%   name, units, and LFP fields
%   Units: units with name and ts fields
%   flags: each flag is names and has ts mat
%   coords: coords for both LED's





% parse varargin
p=inputParser;
addOptional(p,'choose',0);
addOptional(p,'filetype','nex');
parse(p,varargin{:});
choose=p.Results.choose;

%% fetch all NEX files

if ~exist('dirName','var')
    %pick which directory
    dirName=uigetdir('','Pick a directory with one day of nex files');
elseif isempty(dirName)
    dirName=uigetdir('','Pick a directory with one day of nex files');
end





if any(dirName)
    %get all the nex files within
    fileList = getAllNEXFiles(dirName);
    
    dash=WhichDash;
    Session.name=[dirName(find(dirName==dash,1,'last')+1:end) ' Session'];
    %loop through nex and get relevant data
    
    
    % fill out session fields
    %Session.
    % set up some variables
    % a fix just in case you split your tetrode files
    totcel=1;
    % set up an events file

    fprintf('Found the following NEX Files; \n');
    
    for i=1:length(fileList)
        fprintf('%s \n',fileList{i});
    end
else
    fprintf('no nex files found \n');
    return
end


for i=1:length(fileList)
    
    % convert file to struct
    datafile{i} = readNexFileM(fileList{i});
    
    fprintf('Reading %s \n \n', fileList{i}(find(fileList{i}=='\',1,'last')+1:end-4));
    % now go through struct to see whats in it
    
    % get coordinate data if there is any
    if  isfield(datafile{i},'markers');
        % if this file has markers, it may be the tracking coords
        markernames=cellfun(@(a) a.name,datafile{i}.markers,'uni',0);
        % Lets see what we can nab
        if any(cell2mat(strfind(markernames,'DVT')));
            fprintf('found DVT coords \n');
            coords.DVTcoords=readNexCoords(datafile{i});
        elseif any(cell2mat(strfind(markernames,'AVI')));
            fprintf('found AVI coords \n');
            coords.AVIcoords=readNexCoords(datafile{i});
        elseif any(cell2mat(strfind(markernames,'PLX')));
            fprintf('found PLX coords \n');
            coords.PLXcoords=readNexCoords(datafile{i});
        elseif any(cell2mat(strfind(markernames,'Strobe')));
            fprintf('found Strobe coords \n');
            coords.strobe=readNexCoords(datafile{i});
        else
            fprintf('Tracking marker name ''%s'' not recognized \n', markernames{1});
        end
    end
    
    
    % if the file has any LFP data
    if isfield(datafile{i},'contvars')
        fprintf('found LFP data \n');
        Session.LFP.contvars= datafile{i} .contvars;
        Session.LFP.events= datafile{i} .events;
        
    end
    
    
    % if there are neurons, or units
    if isfield(datafile{i},'neurons')
        fprintf('found %d units \n', length(datafile{i}.neurons));
        for j=1:length( datafile{i}.neurons)
            % trim blank space before and after
            tempname= strtrim(datafile{i}.neurons{j}.name);
            % remove all the spaces
            unitdata.units(totcel).name = tempname(tempname~='_') ;
            % get timestamps
            unitdata.units(totcel).ts = datafile{i}.neurons{j}.timestamps;
            % pull in a waveform too
            if isfield(datafile{i},'waves')
                unitdata.units(totcel).wave=datafile{i}.waves{j}.waveforms;
            end
            totcel=totcel+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% interneuron filter goes her %%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    
    % Now go through and see if we have any events we want to grab
    if isfield(datafile{i},'events')
        fprintf('found %d events \n', length(datafile{i}.events));
        
        % tag onto any found events
        rawevents=[rawevents; datafile{i}.events];
    end
end

if ~isfield(Session,'LFP')
    newlfp=input('would you like to import separate LFP data? 1 or 0 \n');
    if newlfp
        dash=WhichDash;
        [lfpfile,lfppath]=uigetfile;
        lfpnex=readNexFileM([lfppath dash lfpfile]);
        if isfield(lfpnex,'contvars')
            Session.LFP.contvars= lfpnex.contvars;
            Session.LFP.events= lfpnex.events;
        else
            Session.LFP=[];
        end
    else
        Session.LFP=[];
    end
end
        % unfortunately the tracking still blows
        Session.coords=coords;
        Session.unitdata=unitdata;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% now process events into design matrices%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('rawevents','var')
    % now make our events
    [flags,pot_matrix,board_matrix,pot_legend,board_legend]=GetDelayEvents(rawevents);
    
    if ~isempty(flags)
        
        % now lets plot out where our samples occurred
        names=fieldnames(coords);
        %keepind=menu({names});
        %for i=1:sum(keepind)
        %end
        tracking=coords.(names{1});
        if ~isempty(tracking)
            sampledata=board_matrix(:,[1 2 5]);
            sampledata(:,3)=sampledata(:,3)+1;
            plotsamples(tracking,sampledata);
        end
        

        Session.flags=flags;
        Session.pot_matrix=pot_matrix;
        Session.board_matrix=board_matrix;
        Session.Board_legend=board_legend;
        Session.pot_legend=pot_legend;
    end
end

% save data out
savedir=uigetdir(dirName,'choose a save folder');
if any(savedir)
    dash=WhichDash;
    save([savedir dash Session.name],'Session');

end

end






function plotsamples(tracking, samples)
figure

% take all timestamps within 200 ms (about 8 frames)
X1=cellfun(@(a) nanmean(tracking(a+.1>(tracking(:,1)) & a+.1< (tracking(:,1)+.05),2)),num2cell(samples(:,1),2));
X2=cellfun(@(a) nanmean(tracking(a+.1>(tracking(:,1)) & a+.1< (tracking(:,1)+.05),4)),num2cell(samples(:,1),2));
Y1=cellfun(@(a) nanmean(tracking(a+.1>(tracking(:,1)) & a+.1< (tracking(:,1)+.05),3)),num2cell(samples(:,1),2));
Y2=cellfun(@(a) nanmean(tracking(a+.1>(tracking(:,1)) & a+.1< (tracking(:,1)+.05),5)),num2cell(samples(:,1),2));
X=nanmean([X1 X2],2);
Y=nanmean([Y1 Y2],2);
pos=[X Y];


% there are four poitions, plot each of them in the positions according to
% the tracking

%
B_R=samples(:,2)==1 & samples(:,3)==1;
plot(pos(B_R,1),pos(B_R,2),'o','color','g')
hold on
B_L=samples(:,2)==2 & samples(:,3)==1;
plot(pos(B_L,1),pos(B_L,2),'x')
P_R=samples(:,2)==1 & samples(:,3)==2;
plot(pos(P_R,1),pos(P_R,2),'x','color','r')
P_L=samples(:,2)==2 & samples(:,3)==2;
plot(pos(P_L,1),pos(P_L,2),'o','color','k')
legend('Beads Right','Beads Left','Rocks Right','Rocks Left');

end



