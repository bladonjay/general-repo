function [Codingstats] = TblocksAnovaFun(samples2,spkmat)
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
    % side
    [p]=anovan(inputparam(:,1),spkmat(:,i),'display','off');
    side(i)=p;
    
    % pot
    [p]=anovan(inputparam(:,2),spkmat(:,i),'display','off');
    pot(i)=p;
    
    % rewarded
    [p]=anovan(inputparam(:,3),spkmat(:,i),'display','off');
    reward(i)=p;
    
    % side by pot, reform into all unique combos
    [~,~,sidepotdes]=unique(inputparam(:,[1 2]),'rows');
    [p]=anovan(sidepotdes,spkmat(:,i),'display','off');
    sidepot(i)=p;
    
    % side by reward
    [~,~,siderewdes]=unique(inputparam(:,[1 3]),'rows');
    [p]=anovan(siderewdes,spkmat(:,i),'display','off');
    siderew(i)=p;
    
    % pot by reward
    [~,~,potrewdes]=unique(inputparam(:,[2 3]),'rows');
    [p]=anovan(potrewdes,spkmat(:,i),'display','off');
    potrew(i)=p;
    
    % first run on rule
    [p]=anovan(rulemat(:,1),spkmat(:,i),'display','off');
    rules(i)=p;
    
    % rule by side
    [~,~,sideruledes]=unique([rulemat(:,1) inputparam(:,1)],'rows');
    [p]=anovan(sideruledes,spkmat(:,i),'display','off');
    siderule(i)=p;
    
    % rule by pot
    [~,~,potruledes]=unique([rulemat(:,1) inputparam(:,2)],'rows');
    [p]=anovan(potruledes,spkmat(:,i),'display','off');
    potrule(i)=p;
    
    % rule by rewarded
    [~,~,rewruledes]=unique([rulemat(:,1) inputparam(:,3)],'rows');
    [p]=anovan(rewruledes,spkmat(:,i),'display','off');
    rewrule(i)=p;
    
    
    
    
    
    % block alone
    [p]=anovan(rulemat(:,2),spkmat(:,i),'display','off');
    % basically pick best, and bonferonni correct
    block(i)=p;
    
    % side by block
    [~,~,sideblockdes]=unique([rulemat(:,2) inputparam(:,1)],'rows');
    [p]=anovan(sideblockdes,spkmat(:,i),'display','off');
    sideblock(i)=p;
    
    % pot by block
    [~,~,potblockdes]=unique([rulemat(:,2) inputparam(:,2)],'rows');
    [p]=anovan(potblockdes,spkmat(:,i),'display','off');
    potblock(i)=p;
    
    % reward by block
    [~,~,rewblockdes]=unique([rulemat(:,2) inputparam(:,3)],'rows');
    [p]=anovan(rewblockdes,spkmat(:,i),'display','off');
    rewblock(i)=p;
    
end
kill







% Struct of all coding types
% first pure coding
Codingstats.rules=rules;
Codingstats.block=block;
Codingstats.side=side;
Codingstats.pot=pot;
Codingstats.reward=reward;

% now conjunctive selectivity
Codingstats.sidepot=sidepot;
Codingstats.siderew=siderew;
Codingstats.potrew=potrew;
Codingstats.potrule=potrule;
Codingstats.siderule=siderule;
Codingstats.rewrule=rewrule;
Codingstats.sideblock=sideblock;
Codingstats.potblock=potblock;
Codingstats.rewblock=rewblock;

end

