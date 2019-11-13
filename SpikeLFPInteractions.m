function [Unitinfo] = SpikeLFPInteractions(Unitinfo, rawLFP, LFPts, FreqBand)
% function [Unit] = SpikeLFPInteractions(Unit, rawLFP, LFPts, FreqBand)
% adds to struct Unit a intrinsic rhythmicity a coherence, and a difference
% in intrinsic vs extrinsic freq power. This can be used to calculate
% precession if it exists.

% I need to split this function into two: autocorrelogram functions and
% spike coherence functions


% backfill all the unit fields:
Unitinfo.burstprob=nan;
Unitinfo.AC.burstindex=nan;
Unitinfo.AC.acthetawelch=nan;
Unitinfo.AC.ACgram=nan(1,141);
Unitinfo.AC.curvestats=nan(2,8);
Unitinfo.AC.thetaRoyer=nan;
Unitinfo.AC.LFPtheta=[nan nan nan];

% i want to do all frequencies... lets see how we can do that
Unitinfo.coherence.coherenceP=nan;
Unitinfo.coherence.coherenceVector=nan;
Unitinfo.coherence.coherenceDeg=nan;



% here are all standard frequency bands to examine
loband=[2 6]; % well there are theta skipping cells and in LEC some low freq cells.
thetaband=[7 12]; % 7-12
betaband=[14 30]; % 25-30 
logammaband=[30 50]; % 40
midgammaband=[50 90]; % 60
higammaband=[90 150]; % 120
rippleband=[150 220];

freqbands=[2 6; 7 12; 14 30; 30 50; 50 90; 90 120; 150 220];
% or just go in 1 or 2 hz increments from 2:220...


% if you dont designate a freq band, theres choices above, but i'll default
% to theta
if ~exist('FreqBand;','var')
    FreqBand=thetaband;
elseif isempty(FreqBand) || isnumeric(FreqBand)
    FreqBand=thetaband;
end


% heres a toggle to nanout units with fewer than 1000 spikes
removebad=0;
% make a plot?

if length(Unitinfo.ts)<100
    warning('This unit has fewer than 100 spikes so stats will be skipped');
end

% if the unit has enough spikes (or if we allow low count units)
if length(Unitinfo.ts)>100 || removebad==1
    % first analyze the autocorrelogram for this frequency
    spikes=Unitinfo.ts;
    timebin=1/1000; % 1 ms per bin
    tseries=accumarray(ceil(spikes/timebin),1); % acgram it at 1 ms steps
    
    % burst index
    % burst index is off of autocorrellogram at 1 millisecond bin lags,
    % and divide the peak between 0 and 10 ms  by the mean value
    % between 40 and 50 ms off
    
    % do a quick acgram for the burst index 
    acgram=xcorr(tseries,tseries,50); 
    % max from 0 to 10 ms lag minus mean 40:50 ms
    % and burst index is ratio of max at short lag to mean at long lag
    if sum(acgram)>100
        Unitinfo.AC.burstindex=max(acgram(40:50))/nanmean(acgram(1:10));
    end

    % burst probability
    
    % a good one is ranck 1973 where they fit the data to a
    % gaussian at -50 to -10 and 10-50 msec and then subtract the
    % guestimate from the real data at -10 to 10 msec
    % andy used 9 ms, royer and buzsaki used 6 ms ISIs
    % burst probability is the % of spikes that have an ISI of less
    % than 9 ms.  So basically get the ISIs and take the % of
    % spikes that have an isi before or after thats <9ms.
    % when i bin i fuck up 60-80 isis, but that really shouldnt
    % impact the overall burst rate
    isi=diff(spikes); % get isis
    closeisi=[isi<.009;0]; % get short pre isis
    closeisi(:,2)=circshift(closeisi,1); % get short post isis
    
    % burst prob is % spikes that have a short pre or post isi
    Unitinfo.burstprob=mean(sum(closeisi,2)>0);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now for fitting the acgram to a parameterized curve
    
    
    % this was taken from a buzsaki paper, the mizuseki
    % paper to analyze theta modulation in dorsal and ventral
    % hippocampus
    
    
    
    % timebin will what fraction of a second?
    timebin=1/100; acgramdist=70;
    step=1000*timebin;
    x=(-step*acgramdist):step:(step*acgramdist);
    tseries=accumarray(ceil(spikes/timebin),1);
    acgram2=xcorr(tseries,tseries,acgramdist);
    % nan out the middle point
    acgram2(acgramdist-1:acgramdist+1)=nan;
    acgram2=acgram2/max(acgram2);
    Unitinfo.AC.ACgram=acgram2;
    % mean out the center spikes (+10, 0, & -10 msec centers) which will bin
    % out to +15 to -15
    acgram2(acgramdist-1:acgramdist+1)=nanmean(acgram2(acgramdist-4:acgramdist+4));

    
    % now model fitting these acgrams to gather theta power
     worsemodel= @(k, ks, kf, tau2, tau3,x)...
        k+... % offset, not sure this is necessary, the slow decay might take care of this
        ks*exp(-abs(x)/tau2)+...              % the slow decay term
        kf*exp(-(x.^2)/tau3^2);              % the fast decay term
    
    % this model is my model with a fast decay a slow decay and a theta
    % add three terms
    bettermodel= @(k, kt, ks, kf, f, tau1, tau2, tau3,x)...
        k+... % offset, not sure this is necessary, the slow decay might take care of this
        kt*cos(.002*pi*f*x).*exp(-abs(x)/tau1)+... % theta with its own decay
        ks*exp(-abs(x)/tau2)+...              % the slow decay term
        kf*exp(-(x.^2)/tau3^2);              % the fast decay term
    % Parameter meanings here:
    % k: overall offset
    % kt: overall theta strength
    % ks: overall strength slow decay
    % kf: overall strength of the fast decay
    % f: frequency of theta or oscillation
    % tau1: decay rate of the oscillation strength
    % tau2: decay rate of the slow bounded by poisson distribution of spiking
    % tau3: decay rate of the bursts
    %  So CF=offset, theta str, burst str, decay str, theta f, theta decay,
    %  slow decay, fast decay
    
    % this is a good curve fitting model, it has four principal
    % bits:
    % ripped from 
    % Distinct representations and theta dynamics in dorsal and ventral hippocampus.
    % from royer et al.
    % 1. a theta term, with a power and a decay rate
    % 2. a slow decay term and offset
    % 3. a fast decay term for bursty cells, or cells twith a very
    % large refractory period (seen in LEC)
    royermodel= @(a, b, c, f, tau1, tau2, x)...
        (a*(cos(.002*pi*f*x)+1)+b).*exp(-abs(x)/tau1)+... % theta with the slow decay
         c*exp(-(x.^2)/tau2^2); % plus fast decay
    
    % parameter meanings here:
    % a: strength of theta
    % b: the slow decay rate less the theta power
    % c: the fast decay strength
    % f: theta frequency
    % tau1: the slow decay rate
    % tau2: the fast decay rate
    
    try % try to fit to this parametric model, it usually works
        
        % really dont need to add my own models here
        %{
        % need to be more intelligent about making htese work, for example
        % the slow decay rates should be bounded to what is expected from a
        % poisson random process
        % the freqbands should also be bounded by the frequency
        %                 [k     ks   kf   tau2   tau3
        [worsemx,wgof] = fit( x', acgram2, worsemodel, ...
            'StartPoint', [.2,   .2,   .2,    100,  50], ...
            'Lower',      [0,   -10,  -10,      10,   0], ...
            'Upper',      [1,    10,   10,    inf, 100],...
            'MaxFunEvals',1000000,'TolFun',10^-10); % 'Robust', 'LAR',
        
         %} 
        [bettermx] = fit( x', acgram2, bettermodel, ...
            'StartPoint', [.2,  1,  .2,   .2,   nanmean(FreqBand),    500,    100,  50], ...
            'Lower',      [0,   0, -10,  -10,   FreqBand(1),    100,      0,   0], ...
            'Upper',      [1,  10,  10,   10,  FreqBand(2),    inf,    inf, 100],...
            'MaxFunEvals',1000000,'TolFun',10^-10); % 'Robust', 'LAR',
       
        %                 [a    b   c    f   tau1   tau2]
        [royermx] = fit(x', acgram2, royermodel, ...
            'StartPoint', [1, .1,   2,   nanmean(FreqBand),  500,  50], ...
            'Lower',      [0,   0, -10,   FreqBand(1),  100,   0], ...
            'Upper',     [10,  10,  10,  FreqBand(2),  inf, 100],...
            'MaxFunEvals',1000000,'TolFun',10^-10); % 'Robust', 'LAR',
        %worsecf=coeffvalues(worsemx);
        bettercf=coeffvalues(bettermx);
        royercf=coeffvalues(royermx);
        % these are basically the three terms
        % see above for information about these
        
        
        
        % this is one way to calculate how much theta power there is in the
        % autocorrelogram
        % sum(abs((x(2)-x(1))*bettercf(2)*exp(-abs(x)/bettercf(6)))))
        
        % the LL is the error basically, and we'll use a normal
        % distribution
        %{
        LL = nansum(log(normpdf(acgram2, modeldata)));
        models(i).AIC = -2*LL + 2*size(mc(nvariables),2);
         % this is just the difference in log likelihood
        models(i).dAIC = models(i).AIC - fullModel.AIC; % difference in AIC
    
        % now grab n params and get the prob of this improving the model
        models(i).n_params = size(ma,2)-size(mc,2); % should just be one...
        models(i).n_params_total = size(mc,2);
        % p val is the chi squared prob of that model over the reduced
        models(i).p = 1-chi2cdf(devC-devA,models(i).n_params);
        %}
        Unitinfo.AC.curvestats(1,1:length(bettercf))=[bettercf];
        Unitinfo.AC.curvestats(2,1:length(royercf))=[royercf];
    catch
      
    end
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Get the spectrogram of the autocorr
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pxx,f] = pwelch(acgram2,50,10,.5:.5:50,1/timebin,'Power');
    %thetapwr=nanmean(pxx(f>5 & f<12));
    
    % normalize by f
    normpower=pxx.*f;
    % upsample to get a good estimate of the peak in the theta range
    upsamplef=.5:.125:50;
    upsamplepwr=interp1(f(f<=50),normpower(f<=50),upsamplef,'spline');
    
    % find peaks
    crossinds=find([0 diff(diff(upsamplepwr)>0)==-1]);
    crossings=[upsamplef(crossinds); upsamplepwr(crossinds)];
    
    % oldschool peaks
    [thetapwr,peakind]=max(upsamplepwr(upsamplef>FreqBand(1) & upsamplef<FreqBand(2)));
    newf=upsamplef(upsamplef>FreqBand(1));
    thetapeak=newf(peakind);

    % get all power above it
    lowpwr=nanmean(normpower(f<FreqBand(1)*0.9));
    hipwr=nanmean(normpower(f>FreqBand(2)*1.1));
    % or all outside the frquency range
    allpwr=normpower(f<FreqBand(1)*0.9 | f>FreqBand(2)*1.1);
    
    % get a z score for how much more theta is above the other power range
    mu=nanmean(allpwr); rho=nanstd(allpwr);
    thetastrength=(thetapwr-mu)/rho;
    
    
    
    Unitinfo.AC.acthetawelch=[thetapwr/nanmean(allpwr) thetapeak];
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% Get the welch spectrum of the AC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if theta power is greater than 1.25% of the other powers in AC
    % then we can calculate its peak with relation to the LFP
    if ~isempty(rawLFP)
    if thetapwr/nanmean(allpwr)>1.25
        % maybe do this for a few epochs?
        [pxx2,f2] = pwelch(rawLFP,1000,500,.5:.5:50,1000,'Power');
        % now grab the peak at theta
        normpower2=pxx2.*f2;
        % now upsample this guy and get its peak
        
        
        upsamplepwr2=interp1(f2(f2<=50),normpower2(f2<=50),upsamplef,'spline');
        
        [thetapwr2,peakind2]=max(upsamplepwr2(upsamplef>FreqBand(1) & upsamplef<FreqBand(2)));
        newf=upsamplef(upsamplef>FreqBand(1));
        thetapeak2=newf(peakind2);
        
        peakdiff=thetapeak-thetapeak2;
        
        Unitinfo.AC.LFPTheta=[thetapeak thetapeak2 peakdiff];
        
    end
    
    
    
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Now do spike phase coherence %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Phase,Amplitude] = GetLFPBand(rawLFP,LFPts,FreqBand);

    spikephase=interp1(LFPts,Phase,spikes,'nearest');
    
    % only analyze when phase is >1std above the mean power
    keepspike=interp1(LFPts,Amplitude,spikes,'nearest');
    
    % take top 75% of amplitude
    cutoff=prctile(Amplitude,25);
    
    % only take timepoints when the amplitude is high
    spikephase(keepspike<cutoff)=[];
    binct=30; % 12 degrees
    
    offset=pi/(binct*2); % ditch edges per bills advice
    okbins=linspace(-pi+offset,pi-offset,binct);
    [spkcounts,b]=histcounts(spikephase,okbins);
    centers=b(1:end-1)+(b(2)-b(1));
    
    % and a rayleigh test
    rstat=circ_rtest(centers,spkcounts);
    % and the vector length
    vlength=circ_r(centers,spkcounts);
    % mean direction
    myangle=circ_mean(centers,spkcounts);
    
    Unitinfo.Theta.coherenceP=rstat;
    Unitinfo.Theta.coherenceVector=vlength;
    Unitinfo.Theta.coherenceDeg=myangle;
    end
    
end

end

