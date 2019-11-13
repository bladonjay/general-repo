function [samples, flags, legend] = ParseTblocksEvents(rawevents)
%function [samples, flags] = ParseTblocksEvents(rawevents)
% this function pulls relevant flags from an imported nex file and parses
% them into a samples design matrix and a flags struct with only relevant
% events in a more specific format.


% for the samples matrix, each row is a different sampling event, and the
% columns are as follows:

% col 1 is the trial number
% col 2 is side (r=1)
% col 3 is pot id (beads=1)
% col 4 is context
% col 5 is pot correct (beads==1)
% col 6 is Block number (1;6) or 1:8
% col 7 is trial number
% col 8 is dig or nodig (dig=1)
% col 9 as the noncontiguous block type (2 for 1 cxt, 4 for 2
% 10th col is the occurrence of that block type (3 or 4 for sess types)
% 11th column was the combo right?
% 12th col is correct or not
% 13th col is the time of the response
% 14th col is a surrogate for performance (arma)
% 15th col is the first sample of that trial
% 16th col is the last sample of that trial
legend={'sample start','side r=1','potID beads=1','context',...
    'which pot is correct','block number','trial num','response dig=1',...
    'block type','block iteration','correct combo','correct response?',...
    'sample end','arma performance','first sample','last sample'};


    % Spatial plots for each epoch
    % 1. get each epoch;
    for i=1:length(rawevents)
        if any(rawevents(i).timestamps)
            eventnames{i}=rawevents(i).name;
            eventtimes{i}=rawevents(i).timestamps;
            fprintf('Found Marker %s \n',eventnames{i});
        end
    end
    
    
    % Tracking data sucks...
    % maybe ask to redo tracking here?
    
    % time to make a giant design matrix
    
    
    % first get all the samples
    % col 1 is the trial number
    % col 2 is side (r=1)
    % col 3 is pot id (beads=1)
    % col 4 is context
    % col 5 is pot correct (1:2)
    % col 6 is Block number (1;6) or 1:8
    % col 7 is trial number
    % col 8 is dig or nodig
    % col 9 as the noncontiguous block type (2 for 1 cxt, 4 for 2
    % 10th col is the occurrence of that block type (3 or 4 for sess types)
    % 11th column was the combo right?
    % 12th col is correct or not
    % 13th col is the time of the response
   
    clear blockflags
    for i=1:length(eventnames)
        
        % block flags will be 1 ts, 2 context, 3 pot rewarded
        if any(strfind(lower(eventnames{i}),'block'))
            blockflags{i}(1,:)=eventtimes{i};
            % if its a block flag (if one context we dont add 1 or 2)
            if any(strfind(lower(eventnames{i}),'2'))
                % in context 1
                blockflags{i}(2,:)=2;
            else
                % in context 2
                blockflags{i}(2,:)=1;
            end
            if any(strfind(lower(eventnames{i}),'beads')) || any(strfind(lower(eventnames{i}),'purple'))
                % beads can be 1 BEADS WINNING IS 1
                blockflags{i}(3,:)=1;
            else
                % or rocks is 2
                blockflags{i}(3,:)=2;
            end
        
        % first get pot samples
        elseif any(strfind(lower(eventnames{i}),'left'))
            % if its the first time weve hit the left side sample
            if any(strfind(lower(eventnames{i}),'beads')) || ...
                    any(strfind(lower(eventnames{i}),'purple'))
                % get media 1
                L_Beads=eventtimes{i};
                % second is side (R==1, L==2)
                L_Beads(:,2)=2;
                % third is the object (beads==1, rocks==2)
                L_Beads(:,3)=1;
            elseif any(strfind(lower(eventnames{i}),'rocks')) || ...
                    any(strfind(lower(eventnames{i}),'straws')) || ...
                    any(strfind(lower(eventnames{i}),'white'))
                % get media 2
                L_Rocks=eventtimes{i};
                L_Rocks(:,2)=2;
                L_Rocks(:,3)=2;
            else
                fprintf('Marker %s not recognized \n', eventnames{i});
            end
            %find right samplse
        elseif any(strfind(lower(eventnames{i}),'right'))
            % code in half of the objects here
            if  any(strfind(lower(eventnames{i}),'beads')) || ...
                    any(strfind(lower(eventnames{i}),'purple'))
                % get name of media 1
                R_Beads = eventtimes{i};
                R_Beads(:,2)=1;
                R_Beads(:,3)=1;
            elseif any(strfind(lower(eventnames{i}),'rocks')) || ...
                    any(strfind(lower(eventnames{i}),'straws')) || ...
                    any(strfind(lower(eventnames{i}),'white'))
                % get media 2
                % get name of  left media 2
                R_Rocks =  eventtimes{i};
                R_Rocks(:,2)=1;
                R_Rocks(:,3)=2;
            else
                fprintf('Marker %s not recognized \n', eventnames{i});
            end
        elseif any(strfind(lower(eventnames{i}),'board down'))
            Board_Down=eventtimes{i};
        elseif any(strfind(lower(eventnames{i}),'board up'))
            Board_Up=eventtimes{i};
        
        elseif any(strfind(lower(eventnames{i}),'dig'))
            if any(strfind(lower(eventnames{i}),'no'))
                no_dig=eventtimes{i};
            else
                dig=eventtimes{i};
            end
        elseif any(strfind(lower(eventnames{i}),'center'))
            center_stem=eventtimes{i};
        elseif any(strfind(lower(eventnames{i}),'frame'))
            frames=eventtimes{i};
        else
            fprintf('Marker %s wasnt used or recognized \n', eventnames{i});
        
        end
    end
    
    %% now analyzing the rest
    
    blockID=sortrows(cell2mat(blockflags)',1);
    % and 4th is the block number in linear order
    blockID(:,4)=1:size(blockID,1);
    
    samples=sortrows([L_Beads;R_Beads;L_Rocks;R_Rocks],1);
    
    fprintf('there were %d samples \n',size(samples,1));
    % mark each sample by its block
    samples(:,4:6)=0;
    for i=1:length(blockID)
        % have to to it like this for some reason
        samples(samples(:,1)>blockID(i,1),4)=blockID(i,2);
        samples(samples(:,1)>blockID(i,1),5)=blockID(i,3);
        samples(samples(:,1)>blockID(i,1),6)=blockID(i,4);
    end
    % 4th col is the context of the block
    % 5th col is which pot is rewarded
    % 6th is the number block over all (6 for short, 8 for long)
    % 7th mark each sample by its trial number (already sorted)
    for i=1:length(Board_Down)
        samples(samples(:,1)>Board_Down(i),7)=i;
    end
    
    % mark 8th col dig or nodig
    dig(:,2)=1; no_dig(:,2)=2;
    end_sample=sortrows([dig;no_dig],1);
    if length(end_sample)~= length(samples(:,1))
        fprintf('samples dont match up to choices! \n');
    end
    samples(:,8)=end_sample(:,2);
    
    % mark 9th col as the noncontiguous block type (2 for 1 cxt, 4 for 2
    % for two contexts;
    % find out how many combos there were (row 4 and 5)
    blockunique=unique(samples(:,[4 5]),'rows');
    % now do block equivalents;
    for i=1:length(blockunique)
        blockidx=(samples(:,[4 5]))==blockunique(i,:);
        samples(sum(blockidx,2)==2,9)=i;
        % now for each of those blocks pull their actual session block
        % 10th col is the occurrence of that block type (3 or 4 for sess
        % types)
        sessblocknum=samples(sum(blockidx,2)==2,6);
        sessunique=unique(sessblocknum);
        for k=1:length(sessunique)
            samples(sum(blockidx,2)==2 & samples(:,6)==sessunique(k),10)=k;
        end
    end
    % now i can pull each block type and find out the iterations
    
    % 11th column was the combo right?
    % 3rd 1=beads, 5th 1=b block
    samples(:,11)=samples(:,3)==samples(:,5);
    
    % 12th col is correct or not
    % 8= 1=dig 2=nodig 11 1=correct 0=incorrect
    samples(:,12)=samples(:,8)~=samples(:,11)+1;
    
    % 13th col is the time of the response
    samples(:,13)=end_sample(:,1);
    
   
    % 14 is does he know or not
    % trying to get a good performance score
    for k=1:length(unique(samples(:,6)))
        thisblock=samples(samples(:,6)==k,:);
        % now rate each trial by performance
        % pull sets of three by performance
        % or pull an arma for the whole block
        % or moving average padded on front
        performance=zeros(size(thisblock,1),1);
        for q=3:size(thisblock,1)
            if sum(thisblock(q-2:q,12))==3
                performance(q:end)=1;
                break
            end
        end
        samples(samples(:,6)==k,14)=performance;
    end
    fprintf('Rat performance: %2.f%% \n',nanmean(samples(:,12))*100);
    % 15 is first sample in trial
    % 16 is last sample in trial
    
    % trying to get a first and last samples
    trialnum=samples(:,7);
    lasts=zeros(length(trialnum),1);
    firsts=zeros(length(trialnum),1);
    for i=1:length(unique(trialnum))
        firstidx=find(trialnum==i,1,'first');
        lastidx=find(trialnum==i,1,'last');
        lasts(lastidx)=1;
        firsts(firstidx)=1;
        
    end
    %
    samples(:,15)=firsts; samples(:,16)=lasts;
    
    flags=struct('blockID',blockID,'board_up',Board_Up,'board_down',...
        Board_Down,'dig',dig,'nodig',no_dig,'beads_right',R_Beads,'beads_left',...
        L_Beads,'straws_right',R_Rocks,'straws_left',L_Rocks);
    if exist('center_stem','var')
        flags.center_stem=center_stem;
    end
    if exist('frames','var')
        flags.frame=frames;
    end
    
    

end

