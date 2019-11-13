function [Codingstats] = TblocksBootstrapFunCtrl(samples2,spkmat)
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
    [p]=BootPval(inputparam(:,1),spkmat(:,i));
    side(i)=p;
    
    % pot
    [p]=BootPval(inputparam(:,2),spkmat(:,i));
    pot(i)=p;
    
    % rewarded
    [p]=BootPval(inputparam(:,3),spkmat(:,i));
    reward(i)=p;
    
    % side by pot, reform into all unique combos
    [~,~,sidepotdes]=unique(inputparam(:,[1 2]),'rows');
    [p]=BootPval(sidepotdes,spkmat(:,i));
    sidepot(i)=p;
    
    
    % first run on rule
    [p]=BootPval(rulemat(:,1),spkmat(:,i));
    rules(i)=p;
    
    % block alone
    [p]=BootPval(rulemat(:,2),spkmat(:,i));
    % basically pick best, and bonferonni correct
    block(i)=p;

    % now we control for either rule or block and see what we get
    
    for r=1:length(unique(rulemat(:,1))) % only two rules
        
        % take cells only for that rule
        inputparamR=inputparam(rulemat(:,1)==r,:);
        spkmatR=spkmat(rulemat(:,1)==r,:);

        % side
        [p]=BootPval(inputparamR(:,1),spkmatR(:,i));
        sideR(i,r)=p;
        
        % pot
        [p]=BootPval(inputparamR(:,2),spkmatR(:,i));
        potR(i,r)=p;
        
        % rewarded
        [p]=BootPval(inputparamR(:,3),spkmatR(:,i));
        rewardR(i,r)=p;
        
        % side by pot, reform into all unique combos
        [~,~,sidepotdes]=unique(inputparamR(:,[1 2]),'rows');
        [p]=BootPval(sidepotdes,spkmatR(:,i));
        sidepotR(i,r)=p;
        
    end
    for b=1:length(unique(rulemat(:,2)))

        % take only trials in taht block
        inputparamB=inputparam(rulemat(:,2)==b,:);
        spkmatB=spkmat(rulemat(:,2)==b,:);

        % side
        [p]=BootPval(inputparamB(:,1),spkmatB(:,i));
        sideB(i,b)=p;
        
        % pot
        [p]=BootPval(inputparamB(:,2),spkmatB(:,i));
        potB(i,b)=p;
        
        % rewarded
        [p]=BootPval(inputparamB(:,3),spkmatB(:,i));
        rewardB(i,b)=p;
        
        % side by pot, reform into all unique combos
        [~,~,sidepotdes]=unique(inputparamB(:,[1 2]),'rows');
        [p]=BootPval(sidepotdes,spkmatB(:,i));
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



function [pval]=BootPval(groups,spikerates)
% so the math is going to be find the group who is the most
% different from all the others, basically i'll run a dprime on the
% spikes of one sample type vs all others, and then bootstrap that

boots=200;
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
    h=histogram(max(abs(effectszb)),30);
    hold on; line([max(abs(effectsz)) max(abs(effectsz))],[0 max(h.Values)],'Color','r');
    if pval>1
        fprintf('asdf');
    end
    kill
end
end





