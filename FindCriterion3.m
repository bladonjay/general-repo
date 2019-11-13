function [first30,middle30,last30,crittrial] = FindCriterion3(session,criterion,trials)
% [first30,middle30,last30,crittrial]=FindCriterion(session,criterion,trials)
%FindCriterion finds the trial at which the rat acheves 'criterion' over 'n trials'
% splits struct into 3 structs, however so far it only carries over the
% trial substruct and the coords fields

% JHB 8-18-14
first30=[];
middle30=[];
last30=[];

correct=[];
for i=1:length(session.trial);
    correct=[correct session.trial(i).correct];
end

for j=16:length(correct)
    percentage=sum(correct(j-trials:j))/trials;
    if percentage>criterion
        crittrial=j;
        break
    end
end
if crittrial>length(correct)-15
    fprintf('rat doesnt reach criterion early enough, crit at %.f \n', crittrial);
elseif crittrial<30
    fprintf('rat reaches criterion too early, crit at %.f \n',crittrial);
else
first30.trial=session.trial(1:30);
last30.trial=session.trial(end-29:end);
middle30.trial=session.trial(crittrial-19:crittrial+10);

first30.coords=session.coords(1:find(session.coords(:,1)<first30.trial(30).start_eat,1,'last'),:);
last30.coords=session.coords(find(session.coords(:,1)<last30.trial(1).enter,1,'last'):end,:);
middle30.coords=session.coords(find(session.coords(:,1)>middle30.trial(1).enter,1,'first')...
    :find(session.coords(:,1)<middle30.trial(end).start_eat,1,'last'),:);

end

end

