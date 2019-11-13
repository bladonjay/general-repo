function [Codingstats] = TblocksAnovaInverse(samples2,spkmat)
% This code appears in Tblocks=SingleUnit, and runs the anovan stat on all
% the comparisons in the Tblocks task

pot=[]; reward=[]; rules=[]; block=[]; side=[];
sidepot=[]; rulepot=[]; rulepos=[]; blockpot=[]; blockpos=[];

% inputparam is side, pot, rewarded
inputparam=samples2(:,[2 3 11]);
% first col is rule
rulemat=samples2(:,5);
% now its block
[~,~,rulemat(:,2)]=unique(samples2(:,6));


% now run stat on each cell
for i=1:size(spkmat,2)

    % across all other parameters is there a code?
    
    % side
    [p]=anovan(spkmat(:,i),inputparam(:,1),'display','off');
    side(i)=p;
    
    % pot
    [p]=anovan(spkmat(:,i),inputparam(:,2),'display','off');
    pot(i)=p;
    
    % rewarded
    [p]=anovan(spkmat(:,i),inputparam(:,3),'display','off');
    reward(i)=p;
    
    % side by pot, reform into all unique combos
    [~,~,sidepotdes]=unique(inputparam(:,[1 2]),'rows');
    [p]=anovan(spkmat(:,i),sidepotdes,'display','off');
    sidepot(i)=p;
    
    
    % first run on rule
    [p]=anovan(spkmat(:,i),rulemat(:,1),'display','off');
    rules(i)=p;
    
    % block alone
    [p]=anovan(spkmat(:,i),rulemat(:,2),'display','off');
    % basically pick best, and bonferonni correct
    block(i)=p;

    % now we control for either rule or block and see what we get
    
    % control for position
    for j=1:length(unique(inputparam(:,1)))
        poscspk=spkmat(inputparam(:,1)==j,:);
        poscrule=rulemat(inputparam(:,1)==j,:);
    
        [p]=anovan(poscspk(:,i),poscrule(:,1),'display','off');
        rulepos(i,j)=p;

        [p]=anovan(poscspk(:,i),poscrule(:,2),'display','off');
        blockpos(i,j)=p;

        % nest object in here too
        potcspk=spkmat(inputparam(:,2)==j,:);
        potcrule=rulemat(inputparam(:,2)==j,:);
    
        [p]=anovan(potcspk(:,i),potcrule(:,1),'display','off');
        rulepot(i,j)=p;

        [p]=anovan(potcspk(:,i),potcrule(:,2),'display','off');
        blockpot(i,j)=p;

    end
    

end
kill

% Struct of all coding types
% first pure coding
Codingstats.rules=rules;
Codingstats.block=block;
Codingstats.side=side;
Codingstats.pot=pot;
Codingstats.reward=reward;
Codingstats.sidepot=sidepot;

% now control for position or pot and bonferonni correct
Codingstats.siderule=2*min(rulepos,[],2)';
Codingstats.potrule=2*min(rulepot,[],2)';
Codingstats.sideblock=2*min(blockpos,[],2)';
Codingstats.potblock=2*min(blockpot,[],2)';

end

