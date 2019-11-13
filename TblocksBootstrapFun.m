function [Codingstats] = TblocksBootstrapFun(samples2,spkmat,nboots)
% This code appears in Tblocks=SingleUnit, and runs the anovan stat on all
% the comparisons in the Tblocks task


if ~exist('nboots','var')
    nboots=10000;
end


inputparam=samples2(:,[2 3 11]);
pot=[]; reward=[]; rules=[]; block=[]; side=[];
sidepot=[]; potrew=[]; siderew=[];
siderule=[]; potrule=[]; rewrule=[];
sideblock=[]; potblock=[]; rewblock=[];


% first col is rule
rulemat=samples2(:,5);
% now its block
[~,~,rulemat(:,2)]=unique(samples2(:,6));
% break out block into single regressors



for i=1:size(spkmat,2)
    % side r vs l
    [p]=BootPval(inputparam(:,1),spkmat(:,i),nboots);
    side(i)=p;
    
    % pot a vs b
    [p]=BootPval(inputparam(:,2),spkmat(:,i),nboots);
    pot(i)=p;
    
    % rewarded or not (rewarded pot)
    [p]=BootPval(inputparam(:,3),spkmat(:,i),nboots);
    reward(i)=p;
    
    % side by pot, reform into all unique combos
    [~,~,sidepotdes]=unique(inputparam(:,[1 2]),'rows');
    [p]=BootPval(sidepotdes,spkmat(:,i),nboots);
    sidepot(i)=p;
    
    % side by reward
    [~,~,siderewdes]=unique(inputparam(:,[1 3]),'rows');
    [p]=BootPval(siderewdes,spkmat(:,i),nboots);
    siderew(i)=p;
    
    % pot by reward
    [~,~,potrewdes]=unique(inputparam(:,[2 3]),'rows');
    [p]=BootPval(potrewdes,spkmat(:,i),nboots);
    potrew(i)=p;
    
    % first run on rule
    [p]=BootPval(rulemat(:,1),spkmat(:,i),nboots);
    rules(i)=p;
    
    % rule by side
    [~,~,sideruledes]=unique([rulemat(:,1) inputparam(:,1)],'rows');
    [p]=BootPval(sideruledes,spkmat(:,i),nboots);
    siderule(i)=p;
     
    % rule by pot
    [~,~,potruledes]=unique([rulemat(:,1) inputparam(:,2)],'rows');
    [p]=BootPval(potruledes,spkmat(:,i),nboots);
    potrule(i)=p;
    
    % rule by rewarded
    [~,~,rewruledes]=unique([rulemat(:,1) inputparam(:,3)],'rows');
    [p]=BootPval(rewruledes,spkmat(:,i),nboots);
    rewrule(i)=p;

    
    
    % block
    [p]=BootPval(rulemat(:,2),spkmat(:,i),nboots);
    % basically pick best, and bonferonni correct
    block(i)=p;
    
    % side by block
    [~,~,sideblockdes]=unique([rulemat(:,2) inputparam(:,1)],'rows');
    [p]=BootPval(sideblockdes,spkmat(:,i),nboots);
    sideblock(i)=p;
    
    % pot by block
    [~,~,potblockdes]=unique([rulemat(:,2) inputparam(:,2)],'rows');
    [p]=BootPval(potblockdes,spkmat(:,i),nboots);
    potblock(i)=p;
    
    % reward by block
    [~,~,rewblockdes]=unique([rulemat(:,2) inputparam(:,3)],'rows');
    [p]=BootPval(rewblockdes,spkmat(:,i),nboots);
    rewblock(i)=p;
    
    
    
    fprintf([num2str(i) '\n']);
end





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

function [pval]=BootPval(groups,spikerates,boots)
% so the math is going to be find the group who is the most
% different from all the others, basically i'll run a dprime on the
% spikes of one sample type vs all others, and then bootstrap that

groupunique=unique(groups);
% get the mean for each group
effectsz=[];
% for each group, find how different that group mean is from all others
for q=1:length(groupunique)
    % max diff in means
    %           mean          spikerates in group        -   mean spikerate         
    %effectsz(q)=nanmean(spikerates(groups==groupunique(q)))-nanmean(spikerates(groups~=groupunique(q)));
    
    effectsz(q)=dprime(spikerates(groups==groupunique(q)),spikerates(groups~=groupunique(q)));
end

% see who has the most different mean, doesnt really matter which,
% because im gonna boot and that will shuffle
effectszb=nan(length(groupunique),boots);
% run n boot permutations of random groupings
parfor boot=1:boots, groupsb(:,boot)=groups(randperm(length(groups))); end

for q=1:length(groupunique)
    parfor boot=1:boots

       % effectszb(q,boot)=nanmean(spikerates(groupsb(:,boot)==groupunique(q)))-nanmean(spikerates(groupsb(:,boot)~=groupunique(q)));
        
        effectszb(q,boot)=dprime(spikerates(groupsb(:,boot)==groupunique(q)),spikerates(groupsb(:,boot)~=groupunique(q)));
    end
end

% how many boots beat that max?
pval=nanmean(max(abs(effectsz))<max(abs(effectszb)));
debug=0;
if any(debug)
    h=histogram(max(abs(effectszb)),30); %#ok<UNRCH>
    hold on; line([max(abs(effectsz)) max(abs(effectsz))],[0 max(h.Values)],'Color','r');
    if pval>1
        fprintf('asdf');
    end
    kill
end
end





