function [Codingstats] = TblocksDprimeFun(samples2,spkmat)
% This code appears in Tblocks=SingleUnit, and runs the dprimedesign stat on all
% the comparisons in the Tblocks task
inputparam=samples2(:,[2 3 11]);
        pot=[]; reward=[]; rules=[]; block=[]; side=[];
        sidepot=[]; potrew=[]; siderew=[]; 
        siderule=[]; potrule=[]; rewrule=[];
        sideblock=[]; potblock=[]; rewblock=[]; 
        
        
        dprimedesign = @(data,design) dprime(data(design==1),data(~design==1));
        
        
        for i=1:size(spkmat,2)
            [p]=dprimedesign(inputparam(:,1),spkmat(:,i));
            side(i)=p;
            
            [p]=dprimedesign(inputparam(:,2),spkmat(:,i));
            pot(i)=p;
            
            [p]=dprimedesign(inputparam(:,3),spkmat(:,i));
            reward(i)=p;
            
            % side by pot, reform into all unique combos
            [~,~,sidepotdes]=unique(inputparam(:,[1 2]),'rows');
            [p]=dprimedesign(sidepotdes,spkmat(:,i));
            sidepot(i)=p;
            
            % side by reward
            [~,~,siderewdes]=unique(inputparam(:,[1 3]),'rows');
            [p]=dprimedesign(siderewdes,spkmat(:,i));
            siderew(i)=p;
            
            % pot by reward
            [~,~,potrewdes]=unique(inputparam(:,[2 3]),'rows');
            [p]=dprimedesign(potrewdes,spkmat(:,i));
            potrew(i)=p;
            
        end
        
        %%%%%%% Block and Rule %%%%%%%%%%%
        % first col is rule
        rulemat=samples2(:,5)==1;
        % now its block
        moreblocks=unique(samples2(:,6));
        % break out block into single regressors
        for i=1:length(moreblocks), rulemat(:,i+1)=samples2(:,6)==moreblocks(i);
        end
        % now run a friedman test on each
        for i=1:size(spkmat,2)
            
            % first run on rule
            [p]=dprimedesign(rulemat(:,1),spkmat(:,i));
            rules(i)=p;
            
            % rule by side
            [~,~,sideruledes]=unique([rulemat(:,1) inputparam(:,1)],'rows');
            [p]=dprimedesign(sideruledes,spkmat(:,i));
            siderule(i)=p;
            
            % rule by rewarded
            [~,~,rewruledes]=unique([rulemat(:,1) inputparam(:,3)],'rows');
            [p]=dprimedesign(rewruledes,spkmat(:,i));
            rewrule(i)=p;
            
            % pot by rule
            [~,~,potruledes]=unique([rulemat(:,1) inputparam(:,2)],'rows');
            [p]=dprimedesign(potruledes,spkmat(:,i));
            potrule(i)=p;
            
            % rule by pot
            
            
            % now run on each regressor
            for j=2:size(rulemat,2)
                [p]=dprimedesign(rulemat(:,j),spkmat(:,i));
                % basically pick best, and bonferonni correct
                block(j-1,i)=p;
                
                % side by block
                [~,~,sideblockdes]=unique([rulemat(:,j) inputparam(:,1)],'rows');
                [p]=dprimedesign(sideblockdes,spkmat(:,i));
                sideblock(j-1,i)=p;
                
                % pot by block
                [~,~,potblockdes]=unique([rulemat(:,j) inputparam(:,2)],'rows');
                [p]=dprimedesign(potblockdes,spkmat(:,i));
                potblock(j-1,i)=p;
                
                % reward by block
                [~,~,rewblockdes]=unique([rulemat(:,j) inputparam(:,3)],'rows');
                [p]=dprimedesign(rewblockdes,spkmat(:,i));
                rewblock(j-1,i)=p;
                
            end 
        kill

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

