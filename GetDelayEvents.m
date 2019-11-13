function [newevents,pot_matrix,board_matrix,pot_legend,board_legend]=GetDelayEvents(events)
% pulls events that I want from the generic events struct that comes from
% the behavioral flags in the NEX file.

% Justification and algorithm
% first, here are the expected flags:

% ball, brush, wood, lego, beads left, beads right, rocks left, rocks
% right, door 1, door 2, board up, treadmill start, treadmill stop, dig,
% nodig, pause, start, stop, pause, return arm

% pause, return arm, treadmill stop, may not be included



% cast events into a cell array;

if ~iscell(events)
    events2=events; clear events;
    for i=1:length(events2)
        events{i}.name=events2(i).name;
        events{i}.timestamps=events2(i).timestamps;
    end
end



%%

newevents=struct; pot_matrix=[]; board_matrix=[]; pot_legend=[]; board_legend=[];

numevs = length(events);

% get the list of event names;
eventnames=[];
for i=1:length(events)
        eventnames{i,1}=events{i}.name;
        eventnames{i,2}=length(events{i}.timestamps);

end
fprintf('Found event names are: \n');
for i=1:length(events)
    fprintf(' \t ''%s'' with %d occurrences \n',eventnames{i,1},eventnames{i,2});
end


for i=1:numevs
    
    %%%%%%%%%%%%%%
    % Get Boards
    %%%%%%%%%%%%%%%
    
    % Output will be Ball, Brush, Wood, Lego
    
    % each board has 3 cols, 1 ts, 2 valence(ball brush==1) and 3 is ID
    % (1:4)
    if any(strfind(lower(eventnames{i}),'ball'))
        % first is timestamp
        ball=events{i}.timestamps;
        % 2nd is valence
        ball(:,2)=1;
        % 3rd is ID
        ball(:,3)=1;
        
    elseif any(strfind(lower(eventnames{i}),'wood')) || ...
            any(strfind(lower(eventnames{i}),'snowflake'))
        wood=events{i}.timestamps;
        wood(:,2)=2;
        wood(:,3)=2;
        
    elseif any(strfind(lower(eventnames{i}),'brush'))
        brush=events{i}.timestamps;
        brush(:,2)=1;
        brush(:,3)=3;
        
    elseif any(strfind(lower(eventnames{i}),'lego'))
        lego=events{i}.timestamps;
        lego(:,2)=2;
        lego(:,3)=4;
        
        
    elseif any(strfind(lower(eventnames{i}),'board up'))
        board_up=events{i}.timestamps;
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% get media for each side
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % output for this will be R_Beads, L_Beads, R_Rocks, L_Rocks
        
    elseif any(strfind(lower(eventnames{i}),'left'))
        % if its the first time weve hit the left side sample
        if any(strfind(lower(eventnames{i}),'bead'))
            
            % get media 1
            L_Beads=events{i}.timestamps;
            % second is side (R==1, L==2)
            L_Beads(:,2)=2;
            % third is the object (beads==1, rocks==2)
            L_Beads(:,3)=1;
            
        elseif any(strfind(lower(eventnames{i}),'rock'))
            
            % get media 2
            L_Rocks=events{i}.timestamps;
            L_Rocks(:,2)=2;
            L_Rocks(:,3)=2;
            
        else
            
            fprintf('Marker %s not recognized \n', eventnames{i});
        end
        
        
        %find right samplse
    elseif any(strfind(lower(eventnames{i}),'right'))
        % code in half of the objects here
        if  any(strfind(lower(eventnames{i}),'bead'))
            
            
            % get name of media 1
            R_Beads = events{i}.timestamps;
            R_Beads(:,2)=1;
            R_Beads(:,3)=1;
            
        elseif any(strfind(lower(eventnames{i}),'rock'))
            
            % get name of  left media 2
            R_Rocks =  events{i}.timestamps;
            R_Rocks(:,2)=1;
            R_Rocks(:,3)=2;
        else
            fprintf('Marker %s not recognized \n', eventnames{i});
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%
        % get outcomes for each side
        %%%%%%%%%%%%%%%%%%%
        
        % no_dig, dig, and return arm (return arm not usually added)
        
        % e.g. end first
    elseif any(strfind(lower(eventnames{i}),'no dig')) ||...
        any(strfind(lower(eventnames{i}),'nodig')) ||...
        any(strfind(lower(eventnames{i}),'end first'))

        no_dig=  events{i}.timestamps;
        no_dig(:,2)=1;
        
        % this is the last choice of the trial
    elseif any(strfind(lower(eventnames{i}),'dig'))

        dig =  events{i}.timestamps;
        dig(:,2)=2;
        % and return arm
    elseif any(strfind(lower(eventnames{i}),'return'))
        return_arm =  events{i}.timestamps;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% get automatic flags %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % outputs: door1, door2, treadmill
        
        
        % if we record the boards
    elseif any(strfind(lower(eventnames{i}),'board 1')) ||...
            any(strfind(lower(eventnames{i}),'door 1'))
        door1 =  events{i}.timestamps;
    
    elseif any(strfind(lower(eventnames{i}),'board 2'))  ||...
            any(strfind(lower(eventnames{i}),'door 2'))
        door2 =  events{i}.timestamps;
        
        % the treadmill start
    elseif any(strfind(lower(eventnames{i}),'treadmill start')) ||...
            any(strfind(lower(eventnames{i}),'keyboard7'))
        treadmill = events{i}.timestamps;
    elseif any(strfind(lower(eventnames{i}),'treadmill stop'))
        treadmillstop = events{i}.timestamps;
        
        % start and stop events
    elseif any(strfind(lower(eventnames{i}),'start'))
        start = events{i}.timestamps;
        
    elseif any(strfind(lower(eventnames{i}),'stop'))
        stop = events{i}.timestamps;
        
    end
    
end

% now lets check if we've filled in all the required variables;
myvars={'ball','brush','wood','lego','board_up','L_Beads','L_Rocks',...
    'R_Beads','R_Rocks','no_dig','dig','treadmill',...
    'start','stop'};
foundvars=ismember(myvars,who);
if sum(foundvars)==length(foundvars)
    fprintf('found all the needed flags \n');
else % if were broken, let em know what went wrong and leave fx
    for i=1:length(foundvars)
        if foundvars(i)==0
            warning('didnt find variable %s ',myvars{i});
            lost=genvarname(myvars{i},who);
            eval([lost '=[]']);
        end
    end
end

%%
%%%%%%%%%%%%% Now add to the events variable %%%%%%%%%%%%
trialconsensus=[length(dig) length(board_up)];
if length(unique(trialconsensus))>1
    fprintf('mismatch in number of dig ( %d ) and board_up ( %d ) \n',trialconsensus(1),trialconsensus(2));
    return
end
if exist('return_arm','var')
    if length(return_arm)~=length(dig)
        warning('bad number of return arm flags');
    end
else
    warning('no return arm flag');
end
% Fill in automatic events
% make sure there arent too many treadmill ts;
if ~isempty(treadmill);
    for i=1:length(dig)
        % they will most likely be the last treadmill before the first sample
        % of each event, because sometimes we start/stop, and sometimes we hit
        % it by accident instead of door 2.
        try
            usetreadmill(i)=treadmill(find(treadmill<dig(i),1,'last'));
        catch
            usetreadmill(i)=treadmill(find(board_up<dig(i),1,'last'));
            warning('GetRelevantEvents:test','Had to use board up flag for treadmill flag %d',i);
        end
    end
    newevents.treadmill=usetreadmill;
end
% door is easy, because our treadmill is 8 seconds
if exist('door1','var'), 
    newevents.door1=door1;
    % sanity check that the door1 variable happens 8s after treadmill
    for i=1:length(usetreadmill)
        doorup(i)=door1(find(door1>usetreadmill(i),1,'first'));
    end
    fprintf('door goes up an average of %d seconds after treadmill \n',mean(doorup-usetreadmill));
elseif exist('treadmill','var');
    newevents.door1=newevents.treadmill+8;
end

if exist('door2','var'), newevents.door2=door2; end
if exist('start','var')
    newevents.start=start; newevents.stop=stop;
end
% fill in board events
newevents.ball=ball;	newevents.wood=wood;  newevents.brush=brush;
newevents.lego=lego;    newevents.board_up=board_up;

% fill in pot events
newevents.L_Beads=L_Beads;   newevents.R_Beads=R_Beads; newevents.no_dig=no_dig;
newevents.L_Rocks=L_Rocks;   newevents.R_Rocks=R_Rocks; newevents.dig=dig;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% now to fill in extra events%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fill in empty events
newevents.trial_correct=[]; newevents.good_choice=[]; newevents.bad_choice=[];
newevents.first_sample=[]; newevents.last_sample=[]; newevents.sample_start=[];
newevents.sample_end=[];  newevents.door1=[];




% now to match choices w answer

% board starts and the ID
% ts, valence(1/2), object (1:4)
boardids=[ball; brush; wood; lego];
% ball and brush are 1, wood and lego are 2
boardids=sortrows(boardids,1);
newevents.board_start=boardids;

% Response ID nodig is 1, dig is 2
% ts, ID(1/2)
responses=[dig; no_dig];
responses=sortrows(responses,1);
newevents.sample_end=responses;

% and pot ID
% 1 is ts, 
% 2 is side(1/2), 
% 3 is pot(1/2)
potids=[R_Beads; R_Rocks; L_Beads; L_Rocks];
potids=sortrows(potids,1);
%%


%%%% for boardmat and samplemat %%%%%%%%
pot_legend={'timestamp','side(r==1)','pot(beads==1)','Trial number',... % 1:4
            'last board ts','last board valence','last board id','board pot match?',... % 5:8
            'rat response(dig=1)','correct action?','when sample stopped',... % 9:12
            'sample number','first sample?','last sample'}; % 13:15
      
% for each SAMPLE
for i=1:length(potids)
    % index the last board (which is trial number)
    lastboard(i)=find(boardids(:,1)<potids(i,1),1,'last');
    % 4th is trial number
    potids(i,4)=lastboard(i);
    % 5th col most recent board ts, 
    % 6th is board valence
    % 7th is the board ID
    potids(i,5:7)=boardids(potids(i,4),1:3);
    % 8th col==1 if its the right pairing, or 0 if wrong pairing
    potids(i,8)=double(potids(i,3)==potids(i,6));
    
    % 9th is what the rat did (1=nodig, 2=dig)
    % get the next choice (1 if nodig, 2 if dig)
    choiceidx(i)= find(responses(:,1)>potids(i,1),1,'first');
    potids(i,9)=responses(choiceidx(i),2)-1;
    % if the choice was the same as matching (0== nomatch and nodig, 1=match
    % and dig)
    newevents.sample_end(i,3)=double(potids(i,8)==potids(i,9));
    % 10th column is if the rat did the right thing (1 if right 0 if wrong(
    potids(i,10)=double(potids(i,8)==potids(i,9));
    % 11th col is when the sample stopped
    potids(i,11)=newevents.sample_end(i,1);

end

% 12th column is the sample number for that trial
% 13 is whether thats the first sample
% 14th is whether its the last sample
[potids(:,12),potids(:,13),potids(:,14)]=SussSamplesPerTrial(potids(:,4));

% fill in
pot_matrix=potids;
newevents.sample_start=potids(:,1:3);

newevents.good_choice=newevents.sample_end(newevents.sample_end(:,3)==1,1);
newevents.bad_choice=newevents.sample_end(newevents.sample_end(:,3)==0,1);

% get first sample.. start, side, pot, end, first, correct 
newevents.first_sample=[newevents.sample_start(newevents.sample_end(:,2)==1,:) ...
    newevents.sample_end(newevents.sample_end(:,2)==1,:)];

% get last samples .. start, side, pot, end, last, correct 
newevents.last_sample=[newevents.sample_start(newevents.sample_end(:,2)==2,:) ...
    newevents.sample_end(newevents.sample_end(:,2)==2,:)];
       
newevents.trial_correct=[newevents.board_start newevents.last_sample(:,3)];

% for each BOARD
% remember:
% first col=board start
% second col=board valene
% 3rd is board id
board_legend={'board timestamp','board valence','board id',...
            'board up','correct?','num samples','last sample end ts'};


for i=1:length(boardids)
    % 4th is board up 
    boardids(i,4)=newevents.board_up(find(newevents.board_up(:,1)>boardids(i,1)...
        ,1,'first'),1);
    % 5th col is whether the trial was correct
    boardids(i,5)=potids(find(potids(:,4)==i,1,'last'),10);
    % 6th col is number of samples for that trial
    boardids(i,6)=length(find(potids(:,4)==i));
    % 7th col is ts when the trial ended
    boardids(i,7)=potids(find(potids(:,4)==i,1,'last'),11);
end
board_matrix=boardids;
    
    

newevents=orderfields(newevents);


end

function [samplesper,first,last]=SussSamplesPerTrial(TrialCounter)
% get it into column
TrialCounter=reshape(TrialCounter,numel(TrialCounter),1);    
first=[1;diff(TrialCounter)==1];
last=[diff(TrialCounter)==1;1];

% first sample==0 middle= 1 last = 2
for i=1:length(TrialCounter)
   samplesper(i)=sum(TrialCounter(1:i)==TrialCounter(i));
end
end
        