function [cor,allcorr,trialstructure,newspkmat]=eventCorrGlobal(spkmat,samples,varargin)
% samples has to have three columns of ones and zeros, and the number of
% rows should equal the number of rows in spkmat.  The samples is basically
% your design matrix from which you categorize your spkmat rates or
% data




% inputs
%   spkmat= the spikemat from sams matrix, it is an n trials by M
%   datapoints(neurons) matrix where each value is the firing rate for that
%   trial (generally speaking from 0 to 1 seconds after stimulus onset)
%   samples= the samples matrix, col 1,2,and 3 should be dichotomous
%   variables, each consisting of 1 vs 0 or 1 vs 2.  Your superordinate
%   variable (in my task it was context) is first, then your other
%   variables follow.  For me, col 1 was context, col 2 was position and
%   col 3 was item

%   'transform'= how do you want to transform your firing rates (zscore?)
%   'cortype'= how do you want to relate your spike rates for trials?  correlate?
%       mahallanobis? cosine distance?
%   'evenodd'= 1 or 0 do you want to only compare even and odd trials, or are
%       even-even and odd-odd comparisons okay?

% OUTPUTS
% cor: is the struct you use to determine your dendrogram and Dprimes
% allcorr: is an organized correlation matrix where each box is the
%   correlation between two trials
% Trialstructure: this is basically a square matrix showing you which
% trials are which
% newspkmat is the transformed spkmat

% Required functions:
%   linearize, calcSim1,
% made by Sam Mackenzie, edited a lot by Jay Bladon

p=inputParser;

% this gives the option to normalize our spike data, most likely we will
% use a x-score or a normalize by max (1 and 3 respectively)
addOptional(p,'transform',1); % default is to zscore
addOptional(p,'inputnames',{'A','B','C'});

% this makes only even and odd comparisions for all the CPI values (should
% be either 0 or 2, 1 just cuts dataset in half.
addOptional(p,'evenOdd',0);

% this is the correlation type, either corr, cosine, mahallanobis or graph,
% see below for the full description
addOptional(p,'cortype',0);  % default is rank correlation

parse(p,varargin{:});
% name your variables
inputnames=p.Results.inputnames;
% transform your spike matrix
transform=p.Results.transform;
% how do you want to compare trials?
cortype=p.Results.cortype;
% do you want to average?
evenodd=p.Results.evenOdd;


if ~exist('cortype','var')
    calctype='corr';
elseif isempty(cortype) || cortype==0
    calctype='corr';
elseif cortype==1
    calctype='mahalanobis';
elseif cortype==2
    calctype='cosine';
elseif cortype==3
    calctype='hamming';
elseif cortype==5
    calctype='Pearson';
end



tr=[];
trialstructure=[];
newspkmat=[];

% this is your samples matrix
% First force each into 1 or 2 because i will be using them as indexes

%%%%%% dont do this, it forces all levels into a binary
% for i=1:3
%     samples(:,i)=double(samples(:,i)==1)+1;
% end

% dont use tr, but its basically the item pos combo
%tr=(samples(:,6)-1)*4+samples(:,8);

%sort by context, position, then item (or param 1,2, and 3)
[samples,b]=sortrows(samples,[1 2 3]);
%[samples,b]=sortrows(samples,[5 6 8]);


% sort the spike matrix the same way
spkmat=spkmat(b,:);

% dont use either of these
avgmat=nan(16,size(spkmat,2));
avgmat1=nan(16,size(spkmat,2));



% log away the original just in case
spkmatold=spkmat; %keep original (ordered like samples)

% zscore to normalize, probably the best option, you zscore so no single
% neuron or datapoint is more important than all the others (large values
% will skew your data, as for correlation you need normally distribitued
% errors.
if transform==1
    spkmat=zscore(spkmat,[],1);
    
    % or do a nonparametric tiedrank
    % this is a decent option
elseif transform==2
    for i=1:length(spkmat(1,:))
        spkmat(:,i)=tiedrank(spkmat(:,i));
    end
    
    % or else normalize by maximum firing rate
    % good option
elseif transform==3
    % for each column
    for i=1:length(spkmat(1,:))
        % each rate as a precent of max
        if max(spkmat(:,i))>0
            spkmat(:,i)=spkmat(:,i)/max(spkmat(:,i));
        end
    end
    
    % log the rate, but watch out for - infinity (when you log 0)
    % probably not a good option
elseif transform==4
    for i=1:length(spkmat(1,:))
        spkmat(:,i)=log(spkmat(:,i));
        idx=isinf(spkmat(:,i));
        spkmat(idx,i)=0;
        clear idx
    end
    % this is only for timeseries data, because it smooths along the first
    % dimension, not the second
elseif transform==5
    for i=1:length(spkmat(:,1))
        spkmat(i,:)=smooth(spkmat(i,:));
    end
    
end


%calc similarity matrix, you can use calcsim to do this as well if you want
%to do mahalanobis or cosine

allcorr=calcSim1(spkmat,calctype);

%allcorr=nancov(spkmat,'pairwise');






% Sam fucks this up and switches position in context... positions 1,3 are
% in context 1 and posiutions 2,4 are in context 2.  This doesnt matter
% except if he decides to hard code corner-context associations
%%%%%% his is for when we force our variables into binary
%allarms = [ 1 1 1; 1 1 2; 1 2 1; 1 2 2; 2 1 1; 2 1 2; 2 2 1; 2 2 2];
% we need to make this adaptable for multiple levels instead of just two
% save out each possible permutation

% this will show you all the available options
allarms=unique(samples,'rows');
cor.allarms=allarms;

% Label the trial types
for i=1:length(allarms(:,1))
    allnames{i,1}=[inputnames{1} '=' num2str(allarms(i,1))...
        ', ' inputnames{2} '=' num2str(allarms(i,2))...
        ', ' inputnames{3} '=' num2str(allarms(i,3))];
end
cor.allnames=allnames;

% now see if we actually have each possible trial type
% for this there ought to be 8 trial types: 2x2x2.
[~,b] = ismember(samples(:,[1 2 3]),allarms,'rows');
existingcombos=unique(b);
existingarms=allarms(existingcombos,:);
existingnames=allnames(existingcombos,:);

cor.existingnames=existingnames;
cor.existingarms=existingarms;



%calculate even and odd assignments, not sure we need it yet though...

%[~,b] = ismember(samples(:,[6 8]),allarms,'rows');
even = [];
odd = [];
for i = 1:length(existingcombos)
    idx = find(b == existingcombos(i));
    if any(idx)
        even = [even idx(1):2:idx(end)];
        odd = [odd idx(1)+1:2:idx(end)];
    end
end

%%
%build average firing for each item and position combination

if evenodd>0
    eo(even) = 1;
    eo(odd) = 2;
    
    
    % for each neuron get average rate for each combo, and split it in half too
    % so you can compare halves
    for j =1:size(spkmat,2)
        
        % for some reason the accumarray fx updated so I cannot build a 4x2 matrix...
        trialave1 = accumarray(samples(:,[1 2 3]),spkmat(:,j),[],@nanmean,nan);
        % therefore this is going to be a 2x2x2 matrix where the first
        % dimension is the context (a) second is position (B) and third is
        % object (c).  We will collapse so its context, position
        tr1 = cat(1,trialave1(:,:,1),trialave1(:,:,2));
        
        % this is because its a 4x2 matrix (and it goes down columns so
        % indexing is transposed)
        
        linearavg= linearize(tr1'); %sort by context, position, item
        cor.spk(:,j) = linearavg(~isnan(linearavg));
        
        trialave2 = accumarray([samples(:,[1 2 3]) eo'],spkmat(:,j),[],@nanmean,nan); %row = place, col = item
        tr2 = cat(1, trialave2(:,:,:,1),trialave2(:,:,:,2));
        linearavg1=linearize(tr2);
        % not sure I want to not take nans, that might break the matrix
        cor.spk1(:,j) = linearavg1(~isnan(linearavg1)); %sort by context, position, item
        
        
    end
    
    % and now to make an even odd baby correlation matrix if you want to
    % average across half your trials of a given type
    
    babycorr=calcSim1(cor.spk1,calctype);
    %imagesc(babycorr);
    % the 45' is the trialaverage compared to itself.
    allcorr=babycorr;
    
    % have to adapt eoarms to the number of levels in each dimension
    eoarms=[[0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1];...
        [0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1];...
        [0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1]];
    eokeep=[];
    for i=1:size(existingarms,1)
        eokeep=[eokeep eoarms(:,existingcombos(i)*2-1) eoarms(:,existingcombos(i)*2)];
    end
    [X1,Y1]=meshgrid(eokeep(1,:)); % first param
    [X2,Y2]=meshgrid(eokeep(2,:));% second
    [X3,Y3]=meshgrid(eokeep(3,:)); % third
    
    trialstructure=((X1==Y1)+(X2==Y2)+(X3==Y3)*.2+checkerboard(1,8));
    
else
    
    
    
    %define comparisons
    [X1,Y1]=meshgrid(allarms(b,1)); %
    [X2,Y2]=meshgrid(allarms(b,2)); %
    [X3,Y3]=meshgrid(allarms(b,3)); %
    
    trialstructure=((X1==Y1)+(X2==Y2)+(X3==Y3)*.2);
    
end

%%
% for all corr
%only keep top half, this is a roundabout way of doing that tho
keep=logical(triu(ones(size(allcorr,1)),1));
allcorr(~keep)=nan;

newspkmat=spkmat;
% this way allows me to
% get the value for each trial, rather than the below that gets a
% linearized matrix that you cant reconstruct

cor.X1=X1; cor.Y1=Y1;
cor.X2=X2; cor.Y2=Y2;
cor.X3=X3; cor.Y3=Y3;
%% X1= item, X2= position X3= Context

%%%%%%%%%%% store averages %%%%%%%%%%%%%

% mean correlation values

% these are the mean values for each param
cor.xCm= nanmean(linearize(allcorr(X3~=Y3))); %
cor.Cm= nanmean(linearize(allcorr(X3==Y3))); %

cor.xAm= nanmean(linearize(allcorr(X1~=Y1))); %
cor.Am= nanmean(linearize(allcorr(X1==Y1))); %

cor.xBm= nanmean(linearize(allcorr(X2~=Y2))); %
cor.Bm= nanmean(linearize(allcorr(X2==Y2))); %



%now two parameters regardless of third
cor.BCm= nanmean(linearize(allcorr(X3==Y3 & X2==Y2))); %
cor.xBCm= nanmean(linearize(allcorr(X3==Y3 & X2~=Y2))); %
cor.BxCm= nanmean(linearize(allcorr(X3~=Y3 & X2==Y2))); %

cor.ACm= nanmean(linearize(allcorr(X3==Y3 & X1==Y1))); %
cor.xACm= nanmean(linearize(allcorr(X3==Y3 & X1~=Y1))); %
cor.AxCm= nanmean(linearize(allcorr(X3~=Y3 & X1==Y1))); %

cor.ABm= nanmean(linearize(allcorr(X2==Y2 & X1==Y1))); %
cor.xABm= nanmean(linearize(allcorr(X2==Y2 & X1~=Y1))); %
cor.AxBm= nanmean(linearize(allcorr(X2~=Y2 & X1==Y1))); %

% now all three
cor.ABCm= nanmean(linearize(allcorr(X3==Y3 & X2==Y2 & X1==Y1))); %

% cross one witin all others
cor.xABCm= nanmean(linearize(allcorr(X3==Y3 & X2==Y2 & X1~=Y1))); %
cor.AxBCm= nanmean(linearize(allcorr(X3==Y3 & X2~=Y2 & X1==Y1))); %
cor.ABxCm= nanmean(linearize(allcorr(X3~=Y3 & X2==Y2 & X1==Y1))); %
% cross two within one
cor.xAxBCm= nanmean(linearize(allcorr(X3==Y3 & X2~=Y2 & X1~=Y1))); %
cor.AxBxCm= nanmean(linearize(allcorr(X3~=Y3 & X2~=Y2 & X1==Y1))); %
cor.xABxCm= nanmean(linearize(allcorr(X3~=Y3 & X2==Y2 & X1~=Y1))); %

% this doesnt really have much sense either but its accounting
cor.xAxBxCm= nanmean(linearize(allcorr(X3~=Y3 & X2~=Y2 & X1~=Y1))); %


%%
%%%%%%%%% each correlation value values %%%%%%%

cor.xA= linearize(allcorr(X1~=Y1)); % DIFF SIDE
cor.A= linearize(allcorr(X1==Y1)); % SAME SIDE

cor.xB= linearize(allcorr(X2~=Y2)); % DIFF NOVELTY
cor.B= linearize(allcorr(X2==Y2)); % SAME NOVELTY

cor.xC= linearize(allcorr(X3~=Y3)); %DIFF SESS
cor.C= linearize(allcorr(X3==Y3)); %SAME SESS



% cant do same position across contexts
%now two parameters across third
cor.BC= linearize(allcorr(X3==Y3 & X2==Y2)); %SAME CON, SAME POS, ANY ITEM
cor.xBC= linearize(allcorr(X3==Y3 & X2~=Y2)); %SAME CON, DIFF POS, ANY ITEM
cor.BxC= linearize(allcorr(X3~=Y3 & X2==Y2)); %SAME CON, DIFF POS, ANY ITEM

cor.AB= linearize(allcorr(X2==Y2 & X1==Y1)); %SAME CON, ANY POS, SAME ITEM
cor.xAB= linearize(allcorr(X2==Y2 & X1~=Y1)); %SAME CON, ANY POS, SAME ITEM
cor.AxB= linearize(allcorr(X2~=Y2 & X1==Y1)); %SAME CON, ANY POS, SAME ITEM

cor.AC= linearize(allcorr(X3==Y3 & X1==Y1)); %SAME CON, ANY POS, SAME ITEM
cor.xAC= linearize(allcorr(X3==Y3 & X1~=Y1)); %SAME CON, ANY POS, SAME ITEM
cor.AxC= linearize(allcorr(X3~=Y3 & X1==Y1)); %SAME CON, ANY POS, SAME ITEM

% now all three
cor.ABC= linearize(allcorr(X3==Y3 & X2==Y2 & X1==Y1)); %SAME CON, SAME POS, SAME ITEM
cor.xABC= linearize(allcorr(X3==Y3 & X2==Y2 & X1~=Y1)); %SAME CON, SAME POS, DIFF ITEM
cor.ABxC= linearize(allcorr(X3~=Y3 & X2==Y2 & X1==Y1)); %DIFF CON, DIFF POS, SAME ITEM
cor.AxBC= linearize(allcorr(X3==Y3 & X2~=Y2 & X1==Y1)); %SAME CON, DIFF POS, SAME ITEM

cor.xAxBC= linearize(allcorr(X3==Y3 & X2~=Y2 & X1~=Y1)); %SAME CON, DIFF POS, DIFF ITEM
cor.AxBxC= linearize(allcorr(X3~=Y3 & X2~=Y2 & X1==Y1)); %DIFF CON, DIFF POS, SAME ITEM
cor.xABxC= linearize(allcorr(X3~=Y3 & X2==Y2 & X1~=Y1)); %DIFF CON, DIFF POS, SAME ITEM

% this doesnt really have much sense either but its accounting
cor.xAxBxC= linearize(allcorr(X3~=Y3 & X2~=Y2 & X1~=Y1)); %DIFF CON, DIFF POS, DIFF ITEM





%% this sets you up to do mahallanobis distance on trial clusters
if cortype==5
    [design]=[1,1,1; 1,1,0; 1,0,1;  1,0,0; 0,1,1; 0,1,0; 0,0,1;  0,0,0;];
    [designreal]=[1,1,1; 1,1,2; 1,3,1;  1,3,2; 2,2,1; 2,2,2; 2,4,1;  2,4,2;];
    for i=1:8 % item pos categories
        % parse out each category one by one
        allcorr{i}=spkmat(samples(:,1)==designreal(i,1) & ...
            samples(:,2)==designreal(i,2) & ...
            samples(:,3)==designreal(i,3));
    end
    %  X1= item, X2=pos, X3=context
    X1=design(:,3); X2=design(:,2); X3=design(:,1);
    trialstructure=design;
end
end


%  samples
% First column is timestamp
% Second column is correct (1) / incorrect (0)
% Third column is left (1) / right (2)
% Fourth column is item set (context pair) (A/B = 1, C/D = 2)
% Fifth column is west (context 1) or east (context 2)
% Sixth column is position (1-4)
% Seventh column is Odor Pair (A/C = 1, B/D = 2)
% Eighth column is Odor (A = 1, B = 2, C = 3, D = 4)
% Ninth column is duration of the sample.
% Tenth column is the day number
% Eleventh column is actual trial #
% Twelth column is the # of samples on that pot
% Thirteenth column is rat being correct (dug on cor. no dig on incorrect)
% Fourteenth column is last sample of that pot for that trial)
% Fifteenth column is whether to exclude the trial
