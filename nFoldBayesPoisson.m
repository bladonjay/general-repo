function [p_correct,fract_correct_real,fract_correct_null,rawstruct] = nFoldBayesPoisson(SpkMat,Classes,nFold,nBoots)
% function [p_correct,fract_correct_real,fract_correct_null] = nFoldBayesPoisson(SpkMat,Classes,nFold,nBoots)
% INPUTS:
%   SpkMat: n trial by M unit matrix of firing rates (any unit of measure
%   is okay, but generally were in Hz I think






% created John H BladonDec 2019
%  Updated 4/28/2021
%%

ClassIDs=unique(Classes);
% ClassIDs=[1 0]; for binary results
cv = cvpartition(length(Classes), 'kfold',nFold);

rawstruct=struct('class',[],'classPriors',[],'realProb',[]);

% xfold crossval
% this will not work if we ahve cells who are silent during the
% odor and silent during the run...
% true decoding
for k=1:nFold
    % so dim 1 is going to be test trial, dim 2 is unit and dim 3
    % is class
    classPriors=nan(1,size(SpkMat,2),length(ClassIDs));
    for j=1:length(ClassIDs)
        classPriors(:,:,j)=mean(SpkMat(Classes==ClassIDs(j) & cv.training(k),:));
    end
    classPriors(classPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
    testmat=SpkMat(cv.test(k),:); % the m trials x n units matrix
    trueOdors=Classes(cv.test(k),:); % M trial ids
    
    %p(x|spikes)=C(prod across units (mean rates)^spikes_i ) * exp( -tau * sum(all rates at x)
    % C  normalizes sum probability to equal 1
    prodmat=prod((classPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
    summat=sum(classPriors,2); % sum across units (mean_i)
    tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
    probmat=prodmat.*exp(-tau.*summat);
    % sum p across potential odors must ==1 (the c term)
    realProb=probmat./sum(probmat,3);
    [~,decoded]=max(squeeze(realProb),[],2);
    p_success(k)=nanmean(trueOdors==ClassIDs(decoded));
    rawstruct(k).classPriors=squeeze(classPriors);
    rawstruct(k).realProb=squeeze(realProb);
    rawstruct(k).trueOdors=trueOdors; % I think this part will work, im not sure
end
fract_correct_real = mean(p_success);

% do shuffle
wb=waitbar(0,'Starting Boots');
for boot = 1:nBoots
    waitbar(boot/nBoots,wb,'Running Boots');
    shuff_odorIDs=Classes(randperm(length(Classes))); % reorder the matrix so it doesnt match odors
    for k=1:nFold
        classPriors=nan(1,size(SpkMat,2),length(ClassIDs));
        for j=1:length(ClassIDs)
            classPriors(:,:,j)=mean(SpkMat(shuff_odorIDs==ClassIDs(j) & cv.training(k),:));
        end
        classPriors(classPriors==0)=realmin; % cast silent cells into really low rates to prevent error out
        testmat=SpkMat(cv.test(k),:); % the m trials x n units matrix
        trueOdors=shuff_odorIDs(cv.test(k),:); % M trial ids
        prodmat=prod((classPriors.^testmat),2); % prod across units i (mean_i)^spikes_i
        summat=sum(classPriors,2); % sum across units (mean_i)
        tau=1; % tau is basically a conversion to n spikes per second (but we've already done that)
        probmat=prodmat.*exp(-tau.*summat);
        % sum p across potential odors must ==1 (the c term)
        realProb=probmat./sum(probmat,3);
        [~,decoded_null]=max(squeeze(realProb),[],2);
        p_success_null(k)=nanmean(trueOdors==ClassIDs(decoded_null));
    end
    fract_correct_null(boot) = mean(p_success_null);
end
 close(wb);
 p_correct=1-normcdf(fract_correct_real,nanmean(fract_correct_null),...
                nanstd(fract_correct_null));
            
end

%{
OLd body from previous code
%}