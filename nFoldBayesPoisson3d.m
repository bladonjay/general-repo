function [predicted_real,predicted_null,realProbMat,nullProbMat] = nFoldBayesPoisson3d(spkmat,nFold,nBoots)
% function [p_correct,fract_correct_real,fract_correct_null] = nFoldBayesPoisson(SpkMat,Classes,nFold,nBoots)
% this function runs an n-fold naive bayesian classifier using a poisson
% distribution.
% it predicts spkmat classes dim 1, given dim 2 units (dimensions) and dim 3 iterations
% this works only if each trial has a sample in each class... although you
% might be able to nan out the classes on each trial that dont occur

% if you dont have all classes each trial, look in jays jadhavcode folder
% under the symanski paper to get the alternative.

% INPUTS:
%   spkmat is m classes x N units x O trials
%   boots is for the null distribution
%   nFold is how many folds (fraction of trials) do you want to test

% OUTPUTS:
%   predicted_real= predictions from real data
%   predicted_null= predictions from null data
%   realProbMat=the posterior probabilities of each option

% this code was borrowed from the code I wrote for symanski et al. 2020

% JHB, 10-16-2020

% crucial input dimensions
Trials=1:size(spkmat,1);
Neurons=1:size(spkmat,2);
Classes=1:size(spkmat,3);

% preallocate outputs
predicted_real=nan(nFold,length(Classes));
predicted_null=nan(nFold,length(Classes),nBoots);
realProbMat=nan(length(Classes),length(Classes),nFold);
nullProbMat=nan(length(Classes),length(Classes),nFold);


% divvy trials into study and test trials for n folds
cv = cvpartition(length(Trials), 'kfold',nFold);
% xfold crossval
% this will not work if we ahve cells who are silent during the
% odor and silent during the run...
% true decoding
decoded=[];
for k=1:nFold
    % so dim 1 is going to be the iteration, dim 2 is unit and dim 3
    % is the class
    
    % first get the priors, which is the mean value across trials fo
    % traiing set
    Priors=permute(nanmean(spkmat(cv.training(k),:,:),1),[3 2 1]); % mean across study trials (Mclasses x Nunits x 1)
    Priors(Priors==0)=realmin; % cast silent cells into really low rates to prevent error out
    testmat=nanmean(spkmat(cv.test(k),:,:),1); % the 1 x N units x Oclasses (real time)
 
    %p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
    % C  normalizes sum probability to equal 1
    
    % prodmat is a M study classes by N test classes
    % so each row is the likelihood of that time given each real time
    % each column then is given that real time, the relative likelihood of
    % each guess
    prodmat=prod((Priors.^testmat),2); % prod across units i (mean_i)^spikes_i
    
    summat=sum(Priors,2); % sum across units (mean_i)
    tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
    probmat=prodmat.*exp(-tau.*summat);
    % sum p across potential odors must ==1 (the c term)
    realProb=probmat./sum(probmat,1);
    [~,decoded(k,:)]=max(squeeze(realProb),[],1); % squeeze dim3 into dim2 (dim2 is now the likelihood of that study class)
    % p success needs to be recalculated to accoutn for dimension shift
    realProbMat(:,:,k)=squeeze(realProb); % the probability of all others
    [~,predicted_real(k,:)]=max(realProbMat(:,:,k));
end



if nBoots>0
    % do shuffle
    wb=waitbar(0,'Starting Boots');
    for boot = 1:nBoots
        waitbar(boot/nBoots,wb,'Running Boots');
        for i=1:size(spkmat,1) % for each trial
            for k=1:size(spkmat,2)% and each neuron
                spkmatShuff(i,k,:)=spkmat(i,k,randperm(length(Classes)));
            end
        end
        for k=1:nFold
            Priors=permute(nanmean(spkmatShuff(cv.training(k),:,:),1),[3 2 1]); % mean across study trials (Mclasses x Nunits x 1)
            Priors(Priors==0)=realmin; % cast silent cells into really low rates to prevent error out
            testmat=nanmean(spkmatShuff(cv.test(k),:,:),1); % the 1 x N units x Oclasses
            prodmat=prod((Priors.^testmat),2); % prod across units i (mean_i)^spikes_i
            summat=sum(Priors,2); % sum across units (mean_i)
            tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
            probmat=prodmat.*exp(-tau.*summat);
            shuffProb=probmat./sum(probmat,1);
            [~,decodedShuff(k,:)]=max(squeeze(shuffProb),[],1); % squeeze dim3 into dim2 (dim2 is now the likelihood of that study class)
            probMatBoot(:,:,k)=squeeze(shuffProb); % the probability of all others
            [~,predicted_null(k,:,boot)]=max(probMatBoot(:,:,k));
        end
        nullProbMat(:,:,boot)=nanmean(probMatBoot,3);
    end
    close(wb);
end
            
end
%}
%{
OLd body from previous code
%}