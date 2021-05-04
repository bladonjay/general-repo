function [p_correct,fract_correct_real,fract_correct_null] = nFoldBayesPoisson2(SpkMat,SpkMat2,Classes,nFold,nBoots)
% so this takes data from two groups of equivalent trial size.  The
% advantage here is that you can cut out trials such that if your two
% spkmats are of contiguous activity patterns, you can avoid using repeated
% datapoints for study and test sets.
%
%
%
%ClassIDs=unique(Classes);
ClassIDs=[1 0];
cv = cvpartition(length(Classes), 'kfold',nFold);
% xfold crossval
% this will not work if we ahve cells who are silent during the
% odor and silent during the run...
% true decoding
for k=1:nFold
    % so dim 1 is going to be test trial, dim 2 is unit and dim 3
    % is right vs left
    odorPriors=mean(SpkMat(Classes==ClassIDs(1) & cv.training(k),:));
    odorPriors(:,:,2)=mean(SpkMat(Classes==ClassIDs(2) & cv.training(k),:)); % 1 by nunits by 2 odors
    odorPriors(odorPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
    testmat=SpkMat(cv.test(k),:); % the m trials x n units matrix
    trueOdors=Classes(cv.test(k),:); % M trial ids
    
    %p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
    % C  normalizes sum probability to equal 1
    prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
    summat=sum(odorPriors,2); % sum across units (mean_i)
    tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
    probmat=prodmat.*exp(-tau.*summat);
    % sum p across potential odors must ==1 (the c term)
    realProb=probmat./sum(probmat,3);
    [~,decoded]=max(squeeze(realProb),[],2);
    p_success(k)=nanmean(trueOdors==ClassIDs(decoded)');
end
fract_correct_real = mean(p_success);

% do shuffle
wb=waitbar(0,'Starting Boots');
for boot = 1:nBoots
    waitbar(boot/nBoots,wb,'Running Boots');
    shuff_odorIDs=Classes(randperm(length(Classes))); % reorder the matrix so it doesnt match odors
    for k=1:nFold
        % use spkmat for study
        odorPriors=mean(SpkMat(shuff_odorIDs==1 & cv.training(k),:));
        odorPriors(:,:,2)=mean(SpkMat(shuff_odorIDs==0 & cv.training(k),:)); % 1 by nunits by 2 odors
        odorPriors(odorPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
        % and spkmat2 for test
        testmat=SpkMat2(cv.test(k),:); % the m trials x n units matrix
        trueOdors=shuff_odorIDs(cv.test(k),:); % M trial ids
        prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
        summat=sum(odorPriors,2); % sum across units (mean_i)
        tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
        probmat=prodmat.*exp(-tau.*summat);
        % sum p across potential odors must ==1 (the c term)
        realProb=probmat./sum(probmat,3);
        [~,decoded_null]=max(squeeze(realProb),[],2);
        p_success_null(k)=nanmean(trueOdors==ClassIDs(decoded_null)');
    end
    fract_correct_null(boot) = mean(p_success_null);
end
 close(wb);
 p_correct=1-normcdf(fract_correct_real,nanmean(fract_correct_null),...
                nanstd(fract_correct_null));
            
end

%{
OLd body from previous code (one sample)

ClassIDs=unique(Classes);
ClassIDs=[1 0];
cv = cvpartition(length(Classes), 'kfold',nFold);
% xfold crossval
% this will not work if we ahve cells who are silent during the
% odor and silent during the run...
% true decoding
for k=1:nFold
    % so dim 1 is going to be test trial, dim 2 is unit and dim 3
    % is right vs left
    odorPriors=mean(SpkMat(Classes==ClassIDs(1) & cv.training(k),:));
    odorPriors(:,:,2)=mean(SpkMat(Classes==ClassIDs(2) & cv.training(k),:)); % 1 by nunits by 2 odors
    odorPriors(odorPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
    testmat=SpkMat(cv.test(k),:); % the m trials x n units matrix
    trueOdors=Classes(cv.test(k),:); % M trial ids
    
    %p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
    % C  normalizes sum probability to equal 1
    prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
    summat=sum(odorPriors,2); % sum across units (mean_i)
    tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
    probmat=prodmat.*exp(-tau.*summat);
    % sum p across potential odors must ==1 (the c term)
    realProb=probmat./sum(probmat,3);
    [~,decoded]=max(squeeze(realProb),[],2);
    p_success(k)=nanmean(trueOdors==ClassIDs(decoded)');
end
fract_correct_real = mean(p_success);

% do shuffle
wb=waitbar(0,'Starting Boots');
for boot = 1:nBoots
    waitbar(boot/nBoots,wb,'Running Boots');
    shuff_odorIDs=Classes(randperm(length(Classes))); % reorder the matrix so it doesnt match odors
    for k=1:nFold
        odorPriors=mean(SpkMat(shuff_odorIDs==1 & cv.training(k),:));
        odorPriors(:,:,2)=mean(SpkMat(shuff_odorIDs==0 & cv.training(k),:)); % 1 by nunits by 2 odors
        odorPriors(odorPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
        testmat=SpkMat(cv.test(k),:); % the m trials x n units matrix
        trueOdors=shuff_odorIDs(cv.test(k),:); % M trial ids
        prodmat=prod((odorPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
        summat=sum(odorPriors,2); % sum across units (mean_i)
        tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
        probmat=prodmat.*exp(-tau.*summat);
        % sum p across potential odors must ==1 (the c term)
        realProb=probmat./sum(probmat,3);
        [~,decoded_null]=max(squeeze(realProb),[],2);
        p_success_null(k)=nanmean(trueOdors==ClassIDs(decoded_null)');
    end
    fract_correct_null(boot) = mean(p_success_null);
end
 close(wb);
 p_correct=1-normcdf(fract_correct_real,nanmean(fract_correct_null),...
                nanstd(fract_correct_null));
            
end

%}