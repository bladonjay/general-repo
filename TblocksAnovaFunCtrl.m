function [Codingstats] = TblocksAnovaFunCtrl(samples2,spkmat)
% This code appears in Tblocks=SingleUnit, and runs the anovan stat on all
% the comparisons in the Tblocks task

pot=[]; reward=[]; rules=[]; block=[]; side=[];
sidepot=[]; potrew=[]; siderew=[];
siderule=[]; potrule=[]; rewrule=[];
sideblock=[]; potblock=[]; rewblock=[];

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
    
    for r=1:length(unique(rulemat(:,1))) % only two rules
        
        inputparamR=inputparam(rulemat(:,1)==r,:);
        spkmatR=spkmat(rulemat(:,1)==r,:);

        % side
        [p]=anovan(spkmatR(:,i),inputparamR(:,1),'display','off');
        sideR(i,r)=p;
        
        % pot
        [p]=anovan(spkmatR(:,i),inputparamR(:,2),'display','off');
        potR(i,r)=p;
        
        % rewarded
        [p]=anovan(spkmatR(:,i),inputparamR(:,3),'display','off');
        rewardR(i,r)=p;
        
        % side by pot, reform into all unique combos
        [~,~,sidepotdes]=unique(inputparamR(:,[1 2]),'rows');
        [p]=anovan(spkmatR(:,i),sidepotdes,'display','off');
        sidepotR(i,r)=p;
        
    end
    for b=1:length(unique(rulemat(:,2)))
        inputparamB=inputparam(rulemat(:,2)==b,:);
        spkmatB=spkmat(rulemat(:,2)==b,:);

        % side
        [p]=anovan(spkmatB(:,i),inputparamB(:,1),'display','off');
        sideB(i,b)=p;
        
        % pot
        [p]=anovan(spkmatB(:,i),inputparamB(:,2),'display','off');
        potB(i,b)=p;
        
        % rewarded
        [p]=anovan(spkmatB(:,i),inputparamB(:,3),'display','off');
        rewardB(i,b)=p;
        
        % side by pot, reform into all unique combos
        [~,~,sidepotdes]=unique(inputparamB(:,[1 2]),'rows');
        [p]=anovan(spkmatB(:,i),sidepotdes,'display','off');
        sidepotB(i,b)=p;
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

% now control for block or rule
Codingstats.sideB=min(sideB,[],2)';
Codingstats.potB=min(potB,[],2)';
Codingstats.rewardB=min(rewardB,[],2)';
Codingstats.sidepotB=min(sidepotB,[],2)';

Codingstats.sideR=min(sideR,[],2)';
Codingstats.potR=min(potR,[],2)';
Codingstats.rewardR=min(rewardR,[],2)';
Codingstats.sidepotR=min(sidepotR,[],2)';

end

