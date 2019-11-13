
clear
%load('C:\Users\Sam\Dropbox\AndreasData\NexFiles-Normal\AJF033\ok\events')

%pick which directory
dirName1=uigetdir('','Pick a directory with one day of nex files');
cii=1;



if any(dirName1)
    
    %get all the nex files within
    fileList = getAllNEXFiles(dirName1);
    
    
    %loop through nex and get relevant data
    totcel=1;
    for i=1:length(fileList)
        nexFile{i} = readNexFile(fileList{i});
        
        
        
        if isfield(nexFile{i},'markers') && isfield(nexFile{i},'events') && (any(cellfun(@(a) strcmp(a.name,'Frame Marker'),nexFile{i}.events)) || any(cellfun(@(a) strcmp(a.name,'Frame_Marker'),nexFile{i}.events)))
            %if the nex file has video coded data save it under cinedata
            cinedata(cii).freq=nexFile{i}.freq;
            cinedata(cii).tbeg=nexFile{i}.tbeg;
            cinedata(cii).tend=nexFile{i}.tend;
            cinedata(cii).events=nexFile{i}.events;
            cinedata(cii).markers=nexFile{i}.markers;
            cinedata(cii).filename=fileList{i};
            
            
            %if it has LFP then save under LFP
            if isfield(nexFile{i},'contvars')
                LFP.contvars= nexFile{i} .contvars;
                LFP.events= nexFile{i} .events;
                
            end
            cii=cii+1;
            
        elseif isfield(nexFile{i},'contvars')  && isfield(nexFile{i},'neurons') && isfield(nexFile{i},'events')
            LFP.contvars= nexFile{i} .contvars;
            LFP.events= nexFile{i} .events;
            %if nex has neuron get neural spike data
            for j=1:length( nexFile{i}.neurons)
                unitdata.units(totcel).name = nexFile{i}.neurons{j}.name;
                unitdata.units(totcel).ts = nexFile{i}.neurons{j}.timestamps;
                totcel=totcel+1;
            end
            
        elseif isfield(nexFile{i},'contvars')
            LFP=nexFile{i};
        elseif isfield(nexFile{i},'neurons')
            for j=1:length( nexFile{i}.neurons)
                unitdata.units(totcel).name = nexFile{i}.neurons{j}.name;
                unitdata.units(totcel).ts = nexFile{i}.neurons{j}.timestamps;
                totcel=totcel+1;
            end
            
        end
    end
end


%%

%get rat info
dash=regexp(dirName1,'\');
options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='tex';

answer = inputdlg({'Rat','date'},dirName1(dash(end)+1:end),1,{'',''},options);

day = answer{2};
rat = answer{1};
%%
% if you ever have two video files or more AND you merged your neural data
% to cut, then this will add the correct time offset into the video markers
% that matches the timeoffset in the neural file


cond=dirName1(dash(end)+1:end);

delta=0;
button = questdlg('Is this session merged?');
%fine the duration of the video coded data based on frame marker

for i=1:length(cinedata)
    frameind{i}=cellfun(@(a) strcmp(a.name,'Frame Marker') | strcmp(a.name,'Frame_Marker'),cinedata(i).events);
    videoLength(i)=cinedata(i).events{ frameind{i}}.timestamps(end);
end
if strcmp(button,'Yes')
    %%realign events to merged spike data
    
    %in the AD channel get start/stop evs
    startind=cellfun(@(a) strcmp(a.name,'Event002'),LFP.events);
    
    if ~any(startind)
        startind=cellfun(@(a) strcmp(a.name,'Start'),LFP.events);
        startEvs=LFP.events{startind}.timestamps;
    else
        startEvs=LFP.events{startind}.timestamps;
    end
    
    stopind=cellfun(@(a) strcmp(a.name,'Stop'),LFP.events);
    
    stopEvs=LFP.events{stopind}.timestamps;
    % dur=stopEvs-startEvs; %duration of each merged session
    
    
    spkLength = diff([startEvs;stopEvs(end)]);
    spkLength(1:end-1)=spkLength(1:end-1)-60;
    [match,indMatch]=bestmatch(videoLength,spkLength);
    
    if any(abs(videoLength-match'))>1
        msgbox('Warning some videos not alligned')
    end
    delta = startEvs(indMatch);
    
    
end

[delta,idx]=sort(delta);
cinedata=cinedata(idx);
frameind=frameind(idx);




frameMarker = [];
for i=1:length(cinedata)
    if any(frameMarker) && cinedata(i-1).events{frameind{i-1}}.timestamps(end)  > (cinedata(i).events{frameind{i}}.timestamps(1) +delta(i))
        frameMarker=[frameMarker;cinedata(i).events{frameind{i}}.timestamps+ cinedata(i-1).events{frameind{i-1}}.timestamps(end)];
    else
        frameMarker=[frameMarker;cinedata(i).events{frameind{i}}.timestamps+delta(i)];
        
    end
end
[~,indDelta]=histc(delta,frameMarker);
%adjust video times by the shift

unitdata.rawLEDs=nan(length(frameMarker),5);
unitdata.rawLEDs(:,1)=frameMarker;
%%

%get Keene video markers labels to make sense
lidx = false;
ridx = false;
numevs = length(cinedata.events);

%loop through all nex events
for i=1:numevs
    
    %find left side samples (one for each pot)
    if any(strfind(lower(cinedata.events{i}.name),'left side'))
        
        if ~lidx
            temp = regexp(cinedata.events{i}.name,'left side');
            medial(1).name = cinedata.events{i}.name(11:end);
            medial(1).ts = cinedata.events{i}.timestamps;
            
            lidx = true;
        else
            
            temp = regexp(cinedata.events{i}.name,'left side');
            medial(2).name = cinedata.events{i}.name(11:end);
            medial(2).ts = cinedata.events{i}.timestamps;
        end
        
        
        %find right samplse
    elseif any(strfind(lower(cinedata.events{i}.name),'right side'))
        
        if ~ridx
            
            
            
            temp = regexp(cinedata.events{i}.name,'right side');
            mediar(1).name = cinedata.events{i}.name(12:end);
            mediar(1).ts = cinedata.events{i}.timestamps;
            
            ridx = true;
            
        else
            temp = regexp(cinedata.events{i}.name,'right side');
            mediar(2).name = cinedata.events{i}.name(12:end);
            mediar(2).ts = cinedata.events{i}.timestamps;
        end
        
    elseif any(strfind(lower(cinedata.events{i}.name),'end of first'))
        endevs(1).name = 'all rejections';
        endevs(1).ts = cinedata.events{i}.timestamps;
        
        
    elseif any(strfind(lower(cinedata.events{i}.name),'incorrect choice'))
        endevs(2).name = 'incorrect';
        endevs(2).ts = cinedata.events{i}.timestamps;
        
        
    elseif any(strfind(lower(cinedata.events{i}.name),'correct rejection'))
        endevs(3).name = 'correct rejection';
        endevs(3).ts = cinedata.events{i}.timestamps;
        
        
    elseif any(strfind(lower(cinedata.events{i}.name),'correct choice'))
        endevs(4).name = 'correct';
        endevs(4).ts = cinedata.events{i}.timestamps;
        
        
        
    elseif any(strfind(lower(cinedata.events{i}.name),'left context enter')) || any(strfind(lower(cinedata.events{i}.name),'west'))
        cinedata.events{i}.name = 'west';
        westts = cinedata.events{i}.timestamps;
    elseif any(strfind(lower(cinedata.events{i}.name),'right context enter')) || any(strfind(lower(cinedata.events{i}.name),'east'))
        cinedata.events{i}.name = 'east';
        eastts = cinedata.events{i}.timestamps;
    elseif any(strfind(lower(cinedata.events{i}.name),'divider up')) || any(strfind(lower(cinedata.events{i}.name),'door'))
        
        cinedata.events{i}.name = 'door';
    elseif any(strfind(lower(cinedata.events{i}.name),'start'))
        
        begints =  cinedata.events{i}.timestamps;
        cinedata.events{i}.name = 'begin';
        
        
        
        
    end
    
end

%%
%match start events with end events
[~,~,c] = unique({medial.name,mediar.name},'stable');
name{1} = medial(1).name;
name{2} = medial(2).name;
%define start identities
%cond = samplesxcondition(column1 = item, column2 = left(1)/right(2)
cond = [...
    c(1)*ones(length(medial(1).ts),1) ones(length(medial(1).ts),1);...
    c(2)*ones(length(medial(2).ts),1)  ones(length(medial(2).ts),1);...
    c(3)*ones(length(mediar(1).ts),1) 2*ones(length(mediar(1).ts),1);...
    c(4)*ones(length(mediar(2).ts),1)  2*ones(length(mediar(2).ts),1)];

tsstart = [medial(1).ts;medial(2).ts;mediar(1).ts;mediar(2).ts];


%sort by ts and maintain condition(cond) label
[tsstart,b] = sort(tsstart);
cond = cond(b,:);

% get all end events
tsend = [ endevs(4).ts ;  endevs(1).ts; endevs(2).ts];

%define end conditions
correct = [ones(length(endevs(4).ts),1); 2*ones(length(endevs(1).ts),1); 3*ones(length(endevs(2).ts),1)];
[tsend,b] = sort(tsend);
correct = correct(b,:);


if length(tsend) ~= length(tsstart)
    disp('End times do not equal start times')
end

if any(diff(tsend - tsstart,[],2)<0)
    disp(['misaligned start/stop time combination ' num2str(round(10*tsstart(diff(tsend - tsstart,[],2)<0))/10)])
end



%divide in to trials
begints = sort(begints);

%figure out if the rat is on east side or west side
[~,trSample] = histc(tsstart,[begints;inf]);
[~,trCon1] = histc(westts,[begints;inf]);
[~,trCon2] = histc(eastts,[begints;inf]);
side = nan(length(trSample),1);
side(ismember(trSample,trCon1)) = 1; %west =1
side(ismember(trSample,trCon2)) = 2; %east =2

%take care of skip trials
while any(isnan(side))
    skip  = find(isnan(side));
    side(skip) = side(skip-1);
end

%figure out which is correct item in correct context
cond = [cond side tsstart];
goodcond = cond(correct==1,:);
[histcond,~,~,b1] = histcn(goodcond(:,[1 3]),1:2,1:2);
[~,b] = max(histcond);


good = [b(1) 1;b(2) 2];
badts = goodcond(~ismember(b1,good,'rows'),:);
for i = 1:size(badts)
    item = name{badts(i,1)};
    
    if badts(i,3)==1
        sid = 'left context';
    else
        sid = 'right context';
    end
    
    disp(['Maybe a mislabeled item? : ' item ' on ' sid ' labeled  correct at ' num2str(round(badts(i,4)*10)/10) 's'])
    
end
%define correct pots for each side
%%
side(side(:,1) == 1,2) = b(1);
side(side(:,1) == 2,2) = b(2);
%%
%define new markers
cinedata.events{numevs+1}.name = 'cl1';
cinedata.events{numevs+1}.timestamps = tsstart(cond(:,2)==1 & cond(:,1) == side(:,2));

cinedata.events{numevs+2}.name = 'il1';
cinedata.events{numevs+2}.timestamps = tsstart(cond(:,2)==1 & cond(:,1) ~= side(:,2));

cinedata.events{numevs+3}.name = 'cr1';
cinedata.events{numevs+3}.timestamps = tsstart(cond(:,2)==2 & cond(:,1) == side(:,2));

cinedata.events{numevs+4}.name = 'ir1';
cinedata.events{numevs+4}.timestamps = tsstart(cond(:,2)==2 & cond(:,1) ~= side(:,2));

cinedata.events{numevs+5}.name = 'cl1end';
cinedata.events{numevs+5}.timestamps = tsend(cond(:,2)==1 & cond(:,1) == side(:,2));

cinedata.events{numevs+6}.name = 'il1end';
cinedata.events{numevs+6}.timestamps = tsend(cond(:,2)==1 & cond(:,1) ~= side(:,2));

cinedata.events{numevs+7}.name = 'cr1end';
cinedata.events{numevs+7}.timestamps = tsend(cond(:,2)==2 & cond(:,1) == side(:,2));

cinedata.events{numevs+8}.name = 'ir1end';
cinedata.events{numevs+8}.timestamps = tsend(cond(:,2)==2 & cond(:,1) ~= side(:,2));



%%



samples=[];
for i=1:length(cinedata)
    temp=multidayeventsKeene(cinedata(i),'before',0);
    temp(:,1)=temp(:,1)+delta(i);
    samples=[samples;temp];
    
    
    
    
    %figure out which position occurs at which time using the frameMarker
    %ts
    [~,bin1]=histc(cinedata(i).markers{1}.timestamps+delta(i),frameMarker);
    [~,bin2]=histc(cinedata(i).markers{2}.timestamps+delta(i),frameMarker);
    
    
    
    
    unitdata.rawLEDs(bin1,2)=cellfun(@(a) str2num(a),cinedata(i).markers{1}.values{1}.strings);
    unitdata.rawLEDs(bin1,3)=cellfun(@(a) str2num(a),cinedata(i).markers{1}.values{2}.strings);
    unitdata.rawLEDs(bin2,4)=cellfun(@(a) str2num(a),cinedata(i).markers{2}.values{1}.strings);
    unitdata.rawLEDs(bin2,5)=cellfun(@(a) str2num(a),cinedata(i).markers{2}.values{2}.strings);
    
    
end

%
% idx = find(diff([0;samples(:,11)])<0);
%
% if any(idx)
%     samples(idx:end,11) = samples(idx:end,11)+samples(idx-1,11);
% end


spkmat=samplespikematrix(samples,unitdata);
%%
%if you have SpksEvs saved already then you can just use this below for the
%SpksEvs to CMB conversion



 dash = regexp(fileList{1},'\');
outfile = [dirName1 '\' rat 'day' num2str(day) '_CMB'];

[t,x,y] = meanCoord(unitdata.rawLEDs);


if size(unitdata.rawLEDs,2)>=5
    hd = headdirectionLED(unitdata.rawLEDs);
    hd(:,2) = hd(:,2).*(180/pi);
    %     [hd, ~, LEDloc] = headdirectionLED(coords);
else
    hd = nan(length(t),2);
end


cinedata.events = cinedata.events(~cellfun(@(a) any(regexp(a.name,'Frame Marker')) |   any(regexp(a.name,'Strobed')),cinedata.events));
event = cell(0,2);
evcount = 0;
for ii = 1:length(cinedata.events)

        newev = numel(cinedata.events{ii}.timestamps);
        event(evcount+1:evcount+newev,1) = {cinedata.events{ii}.name};
        event(evcount+1:evcount+newev,2) = num2cell(cinedata.events{ii}.timestamps);
        evcount = evcount+newev;
    
end



root = CMBHOME.Session('name', outfile, 'b_ts', t,...
    'b_x',x, 'b_y', y, 'b_headdir', hd(:,2),...
    'event', event, 'raw_pos', 1, 'raw_headdir', 1, 'date_created', now, ...
    'epoch', t([1 end]), 'fs_video', 30, 'path_raw_data', fileList{1});
root = root.FixPos;


if isfield(unitdata,'units')
    names = arrayfun(@(x) {x.name}, unitdata.units);
    
    % Get the tetrode numbers for each neuron.
    % This assumes that the tetrode number is the first block of consecutive
    % numbers within the tetrode name.
    tetnum = floor(str2double(regexp(names,'([0-9])+','once','match'))/4)+1;
    
    % Get the unit number for each neuron.
    % This assumes that the units are named such that the last letter of the
    % name denotes the unit number: a = 1, b = 2, c = 3, etc. This is how
    % Offline Sorter names units when using the 'Export to NEX' option.
    unitnum = cellfun(@(x) x(end)-96,names);
    
    for ii = 1:length(unitdata.units)
        fprintf(1,'%s is tetrode %d unit %d\n',names{ii},tetnum(ii),unitnum(ii));
        spk_ts =unitdata.units(ii).ts;
        spk_i = time2ind(t,spk_ts);
        spk_ts = spk_ts(spk_i>0);
        spk_i = spk_i(spk_i>0);
        
        root.spike(tetnum(ii),unitnum(ii)) = CMBHOME.Spike('i',spk_i,'ts',spk_ts,'tet_name',names{ii});
    end
end
%%
%save([dirName1 '\' rat 'day' num2str(day) '_CMB'],'root')



start = [delta indDelta];
figure
X1=cellfun(@(a) nanmean(unitdata.rawLEDs(a+.25>(unitdata.rawLEDs(:,1)) & a+.25< (unitdata.rawLEDs(:,1)+.75),2)),num2cell(samples(:,1),2));
X2=cellfun(@(a) nanmean(unitdata.rawLEDs(a+.25>(unitdata.rawLEDs(:,1)) & a+.25< (unitdata.rawLEDs(:,1)+.75),4)),num2cell(samples(:,1),2));
Y1=cellfun(@(a) nanmean(unitdata.rawLEDs(a+.25>(unitdata.rawLEDs(:,1)) & a+.25< (unitdata.rawLEDs(:,1)+.75),3)),num2cell(samples(:,1),2));
Y2=cellfun(@(a) nanmean(unitdata.rawLEDs(a+.25>(unitdata.rawLEDs(:,1)) & a+.25< (unitdata.rawLEDs(:,1)+.75),5)),num2cell(samples(:,1),2));
X=nanmean([X1 X2],2);
Y=nanmean([Y1 Y2],2);
pos=[X Y];


% there are four poitions, plot each of them in the positions according to
% the tracking
plot(pos(samples(:,6)==4,1),pos(samples(:,6)==4,2),'o','color','g')
hold on
plot(pos(samples(:,6)==1,1),pos(samples(:,6)==1,2),'x')
plot(pos(samples(:,6)==2,1),pos(samples(:,6)==2,2),'x','color','r')
plot(pos(samples(:,6)==3,1),pos(samples(:,6)==3,2),'o','color','k')
CheckPos(X,Y,samples(:,1))
save([dirName1 '\' rat 'day' num2str(day) 'SpksEvs'],'unitdata','spkmat','samples','cinedata','start')
disp(['Saved in: ' dirName1 '\' rat '_' num2str(day) 'SpksEvs.mat'])
