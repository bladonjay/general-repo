function [Ipos]=Positional_Information(X,K)
% function [Ipos]=Positional_Information(X,K);
% uses the positional information metric outlined in Olypher & Fenton 2003
% this metric gathers a positional information score at each bin of k given
% the distribution of firing rates(k) at each bin.

% K and X must be equivalently sized vectors, such that for each X position
% there is a K firing rate.
% generally speaking these are continuous datastreams that are binned into
% ~100 msec bins.
% the math is based on the kullbach lebler divergence
% olypher information scoring
% from olypher and fenton 2003

% the raw spike vectors and matrices
% basic binning parameters;
% 100 Msec time bins
% 5 cm by 5 cm pixels (for a 80 cm arena)

% x and k will be the second and third rows of a huge matrix,

% the input matrix should be somethign like this:
% three column matrix with first column being timebin, second the nspikes,
% and third being the pixel bin.
%** you should filter the bins based on velocity so that you dont capture
%idle time and when the animal is replaying


% JHB 12-7-2020

% xi=pixels: e.g. the linearized positions (will refer the to the unique
% positions that are in the third column of your big matrix
% lets get the occupancy map

Xreal=1:max(X);
for i=1:max(X)
    Xi(i)=sum((X)==i);
end

PXi=Xi/sum(Xi);


% ki= the raw number of spikes happening in each bin
Kreal=0:max(K); % 0 is the min number of spikes
for i=1:max(K)+1
    Ki(i)=sum((K+1)==i); % e.g. at i==1, k=0 which is zero spikes
end

Pki=Ki/sum(Ki);
% Pxi= probability of occupying pixel xi (sum of time bins at each pixel
% over the number of time bins)

Pkxi=Xi.*Ki';
% this is the distribution of firing rates at each pixel
%imagesc(Xreal,Kreal,log(Pkxi)); ylabel('Number of Spikes'); xlabel('Spatial Bin')
% Pk = probability of observing k spikes (how many times did each number of
% spieks occur

% Pxi|k = probability of observing k spikes in pixel xi


% will produce a cell array where each cell is a position bin, and within
% each cell there will be all the time bins (with their rates) that
% occurred in that position

% Pk,xi =  joint probability of observing k spikes and being in pixel xi



for i=1:max(X)+1 % for each spaital bin
    pixelspikes{i}= K(X==(i)); % gather all the timebins there
end

% (e.g. the product of each individual probability)
% and this can actually be a matrix, where the rows are ALL of the observed
% firing rates, and columns are ALL the observed
%sum(ki==(my firing rate) * sum(xi== (my position));

%%  now for the equations
% this will be the vector within each pixel (and you can just take that
% number of spikes divided by all the spikes

% kullbach leibler distances:
% the idea is that this is a distance of one thing to another, its
% directional such that given two distributions (Ai and Bi) across i,
% the distances are:
% d(AB)= sum over i, (A(i)*log(Ai/Bi);


% and then the resistor average is then:
% dAB * dBA / (dAB+dBA) % because these are logs basically

% so to measure positional information across each pixel:

% basically sum over all rates at pixel k, its the prob of that rate at
% that pixel times the log rate at that pixel over all times that rate came
% up

% Ipos(xi) = sum over k Pk|xi log(Pk|xi/Pk)
%          = sum over k Pk|xi (-log(Pk) + log(Pk|xi))

% for each spatial bin
for i=1:length(Xi)
    %tic
    % for each unique spikerate we saw in that bin
    overall=nan;
    myspikes=unique(pixelspikes{i}); % % get all the unique firing rates
    for j=1:length(myspikes)
        % P(k|xi)
        info(j)=sum(pixelspikes{i}==myspikes(j))/length(pixelspikes{i}); % calculate the probability of that spike rate in that bin
        % P(k)
        info2(j)=sum(K==myspikes(j))/length(K); % what % of spikerates overall look like this spike rate
        % log(P(k|xi)/P(k)
        overall(j)=info(j)* log(info(j)/info2(j));
    end
    % sum across all observed k's
    Ipos(i)=nansum(overall);
    %toc
end


% if you mean the abovea cross all pixels, you get the shannons mutual
% information across pixels and position, or I(X,K)


% there are a couple others they calculate:


% additive positional information

% I(a)pos(xi)= - sum across k (Pk * log(Pk)) + sum across k|xi (pk|xi * log(Pk|xi))

% this is nnot too much different except that it can give you negative
% values, and measures the entropy of spike count distributions at each
% pixel vs the entropy at a given pixel


% and finally oyu can calculate the skaggs if you want

%  P=prob at each bin, f_cond is mean # spikes at each bin, and F= mean
%  across bins
%[info, sparsity]=Skaggs_basic(p,f_cond,F);

%%

% so to bin the continuous data:
% first round(timestamps*10) and remove all the zeros
% bincounts=accumarray(round(timestamps*10),1); will get you all the spikes across
% binindices=0:.1:max(round(timestamps,1)) % 100 msec increments
% then you can basically cast the position data into these data by
% interpolating each
%%


% first split the data into bits
plotit=0;
if plotit
    figure; plot(cellfun(@(a) mean(a), pixelspikes)/nanmean(Xi));
    xlabel('Binned position'); ylabel('firing rate');
    yyaxis right;
    plot(Ipos); ylabel('Positional INformation');
end

end



