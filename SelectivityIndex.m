function [score,pval] = SelectivityIndex(rates,design,numboots)
% Score = SelectivityIndex(rates,design)
% maximum selectivity index for a level to a variable over all others
% basically ARGMAX i->n ( abs (a(i)-a(~i)) / (a(i)+a(~i)) )
design=design(:); rates=rates(:); pval=nan; score=nan;
if ~exist('numboots','var')
    numboots=0;
elseif isempty(numboots)
    numboots=0;
end
% for each level, subtract its rate from other levels
score=SIquick(design,rates);

if numboots>0
    bootscore=nan(numboots,1);
    % preallocate the shuffles
    for i=1:numboots, shuffled(i,:)=randperm(length(design)); end
    parfor boot=1:numboots
       bootdesign{boot}=design(shuffled(boot,:));
       bootscore(boot)=SIquick(bootdesign{boot},rates);
    end
    % now the pval is just how many boots beat or tie the real score
    %pval=nanmean(bootscore>=score);
    % other way is to get a normpdf (all vals are abs already)
    if nanstd(bootscore)>0
        pval = 1-normcdf( score, nanmean(bootscore), nanstd(bootscore));
    else
        pval=nan;
    end
    if pval<.05
        fprintf('got one \n');
    end
end
end

function score=SIquick(design,rates)
       [uniques,~,ia]=unique(design);
       scores=nan(length(uniques),1);
       % run for each level of that factor
       for i=1:length(scores)
           % absolute val of difference
           numer=abs(nanmean(rates(ia==i))-nanmean(rates(ia~=i)));
           % sum
           denomer=nanmean(rates(ia==i))+nanmean(rates(ia~=i));
           % absolute difference over the sum
           scores(i)=numer/denomer; 
       end
       % get argmax of those scores
       score=max(scores);
end


