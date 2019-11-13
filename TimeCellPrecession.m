function [precessionstats,precessionimage] = TimeCellPrecession(trspikes,phases)
% function [precessionstats,precessionimage] = TimeCellPrecession(spikes,phases,bins)
% spikes=vector of event locked timestamps
% phases=the phase of those spikes
% bins (not necessary)
% ouitput
%% this searches for fields and finds the field size that gives you the
% most reliable slope for your precession.
% requires corrC2Lin_Kempter2012 and maximizes the -log2(P)
precessionimage=[];
precessionstats.p=nan;
precessionstats.slope=nan;
precessionstats.phasestart=nan;
precessionstats.onedge=1;
precessionstats.fieldwidth=nan;


if length(phases)>50
    
    delayspikes=cell2mat(trspikes');
    meancurve=SmoothMat2(histcounts(delayspikes,0:.001:8),[1000 1000],100);
    % now get peak
    [fieldctr,ctrpos]=max(meancurve); % in msecs
    
    
    % now pan throuigh top 5, 10, 15, 20 % of firing rates to grab fields
    thresh=[95 90 85 80 75 60];
    verbose=0;
    
    % first do the whole delay
    [precessionstats.slope,precessionstats.phasestart,~,p]=...
        corrC2Lin_Kempter2012(delayspikes,phases);
    precessionstats.p=-log2(p);
    precessionstats.fieldwidth(1)=max(delayspikes)-min(delayspikes); % whole delay width
    
    % now progressively step down from the peak yes we're only taking the
    % best peak
    for i=1:length(thresh)
        % last value that was below thresh before the real peak
        thisstart=find(meancurve(1:ctrpos)<prctile(meancurve,thresh(i)),1,'last');
        if isempty(thisstart), thisstart=1; onedge(i+1)=1; end
        % first value after peak
        thisend=ctrpos+find(meancurve(ctrpos:end)<prctile(meancurve,thresh(i)),1,'first');
        if isempty(thisend), thisend=length(meancurve); onedge(i+1)=1;  end
        % now corr the values there
        infield=delayspikes>thisstart/1000 & delayspikes<thisend/1000;
        [precessionstats.slope(i+1),~,~,thisp,ev]=...
            corrC2Lin_Kempter2012(delayspikes(infield)-thisstart/1000,phases(infield));
        precessionstats.p(i+1)=-log2(thisp);

        % grab the expected value of the first spike
        [~,firstspike]=min(delayspikes(infield));
        precessionstats.phasestart(i+1)=ev(firstspike);
        precessionstats.fieldwidth(i+1)=(thisend-thisstart)/1000; % seconds wide
        if verbose && precessionstats.p(i+1)>.98 && precessionstats.p(i+1)<1
            figure; plot(meancurve);
            hold on;
            plot([thisstart thisstart],[0 meancurve(thisstart)],'r');
            plot([thisend thisend],[0 meancurve(thisend)],'r')
            yyaxis right;
            plot(delayspikes*1000,phases,'.'); hold on;
            plot(delayspikes(infield)*1000,ev,'b','LineWidth',2);
            title(sprintf(' p=%.2f slope=%.2f',precessionstats.p(i+1),precessionstats.slope(i+1)));
        end
    end
end

% from there you can find the slope that yielded the lowest p value
end
