function samples = andreaeventsKeene(cinedata, daynum, before, after)
% this organizes the cinedata from the nex struct into a samples matrix
% daynum not a big deal here, but cinedata should have 




if(exist('daynum','var')~=1); daynum = 1; end
if(exist('before','var')~=1); before = []; end
if(exist('after','var')~=1); after = []; end

%cinedata is a struct from readnexfile produced form cineplex
%daynum is the day of the taining session


% First column is timestamp
% Second column is correct (1) / incorrect (0) (pot)
% Third column is left (1) / right (2)
% Fourth column is item set (context pair) (A/B = 1, C/D = 2)
% Fifth column is west (context 1) or east (context 2)
% Sixth column is position (1-4)
% Seventh column is Odor Pair (A/C = 1, B/D = 2)
% Eighth column is Odor (A = 1, B = 2, C = 3, D = 4)
% Ninth column is duration of the sample.
% Tenth column is the day number
% Eleventh column is actual trial #
% Twelth column is the # of samples on that pot
% Thirteenth column is rat being correct (dug on cor. no dig on incorrect)
% Fourteenth column is last sample of that pot for that trial)
% Fifteenth column is whether to exclude the trial

% Initialize sample matrices as empty.
correctleft1 = zeros(0,1);
correctright1 = zeros(0,1);
incorrectleft1 = zeros(0,1);
incorrectright1 = zeros(0,1);
correctleft2 = zeros(0,1);
correctright2 = zeros(0,1);
incorrectleft2 = zeros(0,1);
incorrectright2 = zeros(0,1);

% Get the end time of each sampling
correctleft1end = zeros(0,1);
correctright1end = zeros(0,1);
incorrectleft1end = zeros(0,1);
incorrectright1end = zeros(0,1);
correctleft2end = zeros(0,1);
correctright2end = zeros(0,1);
incorrectleft2end = zeros(0,1);
incorrectright2end = zeros(0,1);


%% Find what event number is used for each event type in the NEX file cell array.
for ii = 1:numel(cinedata.events)
    switch lower(cinedata.events{ii}.name)
        % east and west are pretty consistent, those always pass thru
        case 'west'; west = cinedata.events{ii}.timestamps;
        case 'east'; east = cinedata.events{ii}.timestamps;
        % these are all coded in from mergeunitwithvideokeene
        case 'cl1'; correctleft1 = cinedata.events{ii}.timestamps;
        case 'cr1'; correctright1 = cinedata.events{ii}.timestamps;
        case 'il1'; incorrectleft1 = cinedata.events{ii}.timestamps;
        case 'ir1'; incorrectright1 = cinedata.events{ii}.timestamps;
        case 'cl2'; correctleft2 = cinedata.events{ii}.timestamps;
        case 'cr2'; correctright2 = cinedata.events{ii}.timestamps;
        case 'il2'; incorrectleft2 = cinedata.events{ii}.timestamps;
        case 'ir2'; incorrectright2 = cinedata.events{ii}.timestamps;
        case 'cl1end'; correctleft1end = cinedata.events{ii}.timestamps;
        case 'cr1end'; correctright1end = cinedata.events{ii}.timestamps;
        case 'il1end'; incorrectleft1end = cinedata.events{ii}.timestamps;
        case 'ir1end'; incorrectright1end = cinedata.events{ii}.timestamps;
        case 'cl2end'; correctleft2end = cinedata.events{ii}.timestamps;
        case 'cr2end'; correctright2end = cinedata.events{ii}.timestamps;
        case 'il2end'; incorrectleft2end = cinedata.events{ii}.timestamps;
        case 'ir2end'; incorrectright2end = cinedata.events{ii}.timestamps;
        % these dont pass through, there are none
        case 'door'; door = cinedata.events{ii}.timestamps;
        case 'begin'; begin = cinedata.events{ii}.timestamps;
    end
end

%% Label 'west' (context 1) and 'east' (context 2) events then combine the lists.
west(:,2) = 1;
east(:,2) = 2;
context = sortrows([west; east],1);

% Label sample events
% Second column = correct (1) / incorrect (0)
correctleft1(:,2) = 1;
correctright1(:,2) = 1;
incorrectleft1(:,2) = 0;
incorrectright1(:,2) = 0;
correctleft2(:,2) = 1;
correctright2(:,2) = 1;
incorrectleft2(:,2) = 0;
incorrectright2(:,2) = 0;

% Third column = left (1) / right (2)
correctleft1(:,3) = 1;
correctright1(:,3) = 2;
incorrectleft1(:,3) = 1;
incorrectright1(:,3) = 2;
correctleft2(:,3) = 1;
correctright2(:,3) = 2;
incorrectleft2(:,3) = 1;
incorrectright2(:,3) = 2;

% Fourth column is item set (context pair) (A/B = 1, C/D = 2)
correctleft1(:,4) = 1;
correctright1(:,4) = 1;
incorrectleft1(:,4) = 1;
incorrectright1(:,4) = 1;
correctleft2(:,4) = 2;
correctright2(:,4) = 2;
incorrectleft2(:,4) = 2;
incorrectright2(:,4) = 2;

% Combine the event lists
[samples, order] = sortrows([correctleft1;
    correctright1;
    incorrectleft1;
    incorrectright1;
    correctleft2;
    correctright2;
    incorrectleft2;
    incorrectright2],1);
samplesend = [correctleft1end;
    correctright1end;
    incorrectleft1end;
    incorrectright1end;
    correctleft2end;
    correctright2end;
    incorrectleft2end;
    incorrectright2end];
if any(samplesend)
    samplesend = samplesend(order,:);
else
    samplesend=nan(length(order),1); % not videocoded
    after=1.5; %fiexed time
end

if(numel(order)~=size(samplesend,1))
    error('Number of start timestamps does not match number of stop timestamps.');
    % elseif(~isequal(order,orderend))
    %     error('Order of the stop timestamps does not match the order of the start timestamps.');
elseif(any(samplesend(:,1)<samples(:,1)))
    ind = find(samplesend(:,1)<samples(:,1),1);
    fprintf(2,'Correct/Incorrect: %d, Left/Right: %d, Odor Pair: %d\n',samples(ind,2:4));
    warning('andreaevents:tswarning','Sample end before samples start at timestamp %d (%d).',samples(ind,1),ind);
elseif(any(samplesend(1:end-1,1)>samples(2:end,1)))
    ind = find(samplesend(1:end-1,1)>samples(2:end,1),1);
    fprintf(2,'Correct/Incorrect: %d, Left/Right: %d, Odor Pair: %d\n',samples(ind,2:4));
    warning('andreaevents:tswarning','Sample end overlaps with start of next sample at timestamp %d (%d).',samplesend(ind,1),ind);
end

% Determine in what context each sample occurred.
% Fifth column is west (context 1) or east (context 2)
for ii = 1:size(context,1)-1
    thiscontext = (samples(:,1)>=context(ii,1) & samples(:,1)<context(ii+1,1));
    samples(thiscontext,5) = context(ii,2);
end
thiscontext = (samples(:,1)>=context(end,1));
samples(thiscontext,5) = context(end,2);

%  WEST (1)      EAST (2)
% +--------+    +--------+
% | Pos 1  |    |  Pos 2 |
% |        +----+        |
% |                      |
% |        +----+        |
% | Pos 3  |    |  Pos 4 |
% +--------+    +--------+

% Determine in which position each sample occured.
% Third column is left (1) / right (2)
% Fifth column is west (context 1) or east (context 2)

% Position 1 = Context 1 on the Right (2)
pos1 = (samples(:,5)==1 & samples(:,3)==2);
samples(pos1,6) = 1;

% Position 2 = Context 2 on the Left (1)
pos2 = (samples(:,5)==2 & samples(:,3)==1);
samples(pos2,6) = 2;

% Position 3 = Context 1 on the Left (1)
pos3 = (samples(:,5)==1 & samples(:,3)==1);
samples(pos3,6) = 3;

% Position 4 = Context 2 on the Right (2)
pos4 = (samples(:,5)==2 & samples(:,3)==2);
samples(pos4,6) = 4;

% Determine which odor was sampled
% Second column is correct (1) / incorrect (0)
% Fifth column is west (context 1) or east (context 2)
% Seventh column is Odor Pair (A/C = 1, B/D = 2)

% Odor Pair 1 (A or C) was correct in west (context 1)
item1 = (samples(:,5)==1 & samples(:,2) == 1);
samples(item1,7) = 1;

% Odor Pair 2 (B or D) was incorrect in west (context 1)
item2 = (samples(:,5)==1 & samples(:,2) == 0);
samples(item2,7) = 2;

% Odor Pair 1 (A or C) was incorrect in east (context 2)
item1 = (samples(:,5)==2 & samples(:,2) == 0);
samples(item1,7) = 1;

% Odor Pair 2 (B or D) was correct in east (context 2)
item2 = (samples(:,5)==2 & samples(:,2) == 1);
samples(item2,7) = 2;

% Identify the individual odors
% Fourth column is item set (context pair) (A/B = 1, C/D = 2)
% Seventh column is Odor Pair (A/C = 1, B/D = 2)
% Eighth column is Odor (A = 1, B = 2, C = 3, D = 4)

% Odor A occurs in Set 1 and Odor Pair 1
odorA = (samples(:,4)==1 & samples(:,7) == 1);
samples(odorA,8) = 1;

% Odor B occurs in Set 1 and Odor Pair 2
odorB = (samples(:,4)==1 & samples(:,7) == 2);
samples(odorB,8) = 2;

% Odor C occurs in Set 2 and Odor Pair 1
odorC = (samples(:,4)==2 & samples(:,7) == 1);
samples(odorC,8) = 3;

% Odor D occurs in Set 2 and Odor Pair 2
odorD = (samples(:,4)==2 & samples(:,7) == 2);
samples(odorD,8) = 4;

% Ninth column of samples is the duration of the sample: default actual
% duration
samples(:,9) = samplesend(:,1)-samples(:,1);

% Tenth column is the day number
samples(:,10) = daynum;


%override the default with fixed window length
if(~isempty(after) && ~isempty(before))
    if(before == 0 && after == 0)
        error('andreaevents:ZeroDuration','Before and after are both zero, meaning zero duration.');
    end
    warning('andreaevents:BeforeandAfter',...
        'Using both %d seconds before and %d seconds after, ignoring the sample end times.',before,after);
    samples(:,1) = samples(:,1)-before;
    samples(:,9) = before+after;
elseif(~isempty(before));
    warning('andreaevents:OnlyBefore',...
        'Using %d seconds before, but keeping sample end times.',before);
    samples(:,1) = samples(:,1)-before;
    samples(:,9) = samples(:,9)+before;
elseif(~isempty(after));
    warning('andreaevents:OnlyAfter',...
        'Using %d seconds after sample end, but keeping sample start times.',after);
    samples(:,9) = samples(:,9)+after;
end
if any(samples)
    samples = fixSamples(samples, door,begin,context);
end
end
