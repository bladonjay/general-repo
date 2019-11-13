function [newts,oldts,corvals,offsettime]=AlignBySpikes(ts,unitdata,style)

% The idea here is that you can get better time resolution of your events
% if you use the spike trains to inform the events.  The idea is that the
% temporal structure of the aggregate tuning curve is more accurate than
% our video coding.  Therefore we can adjust our video coding a little bit
% based on the spike train.

% INPUTS
% ts: vector of timestamps corresponding to estimated event times
% unitdata: struct with units field each units field has 'ts'
% style: 1 or 2, see below for basic algorithm
% Basic Algorithm:
%   1. For Sequential Alignment:
%       a. build a tuning curve of each cell centered around the event for
%       trial 1 and trial 2 and smooth (probably based on the average
%       change between sequential ts)
%       b. adjust the ts for event 2 to the max correlation between events
%       c. do the same for the next trial
%   2. For recursive to Average:
%       a. build an AVERAGE tuning curve for each cell centered around the
%       event.
%       b. Find max correlation of each trial aggregate tuning curve to the
%       average.
%       c. optional: recalculate average after each trial is adjusted
% OUTPUTS
% newts: the new timestamps that are lined up based on spike data
% oldts: the old ts, so you can compare
% corvals: these are the correlation values of each event to the average
% offsettime: these are the offsets themselves
%%

% For first method:
oldts=ts;
newts=oldts;
corvals=[]; offsettime=[];

% first develop the timebins we're gonna use, lets go from -.4 to 4s in 4
% ms bins, that way we smooth over about a 20 ms window, enough to get
% within a quarter of a theta cycle and hopefully grab gamma peaks
binsz=.002;
bins=-.1:binsz:.1;

if style==1
    h=waitbar(0,'lining up the trials');
    for i=2:length(ts)
        % get the timebins for this event
        trialbints1=bins+newts(i-1);
        % get the tuning curve for the first event
        smooth1s=makecurve(trialbints1,unitdata,binsz);
        
        % get the timebins for second event
        trialbints2=bins+newts(i);
        % get the tuning curve for the second event
        smooth2s=makecurve(trialbints2,unitdata,binsz);
        
        
        % gets you a 2d cross cor (using a convolution)
        % and we can eliminate edge effects by nan'ing the corrs where there are
        % only a few overlapping pixels, say 5x the number of bins there are
        [acmat]=xcorr2(smooth1s,smooth2s);
        
        % the ac will double the size of our bins becasuse it will convolve
        actimebins=min(bins)*2:binsz:max(bins)*2;
        % we want a zero x offset and to find the max at the y offset
        [~,cols]=size(acmat);
        centercol=ceil(cols/2);
        % take only the center column (where the cells line up) and find the max.
        % the max can be far away from the center, so you may want to convolve the
        % plot with a cosine?
        offsetcors=acmat(:,centercol);
        
        [maxcorr,offsetind]=max(offsetcors);
        offsettime(i)=actimebins(offsetind);
        newts(i)=newts(i)-offsettime(i);
        corvals(i)=maxcorr;
        waitbar(i/length(ts))
    end
    close(h)
end
%%
% for second method
if style==2
    fprintf('making average tuning curves \n ');
    for i=1:length(ts)
        allcurves(:,:,i)=makecurve(ts(i)+bins,unitdata,binsz);
    end
    fprintf('done making average curves \n');
    
    h=waitbar(0,'lining up the trials');
    for i=1:length(ts)
        
        meancurves=mean(allcurves,3);
        
        % get the timebins for this event
        trialbints1=bins+newts(i);
        % get the tuning curve for the first event
        smooth1s=makecurve(trialbints1,unitdata,binsz);
        
        [acmat]=xcorr2(meancurves,smooth1s);
        
        % the ac will double the size of our bins becasuse it will convolve
        actimebins=min(bins)*2:binsz:max(bins)*2;
        % we want a zero x offset and to find the max at the y offset
        [~,cols]=size(acmat);
        centercol=ceil(cols/2);
        % take only the center column (where the cells line up) and find the max.
        % the max can be far away from the center, so you may want to convolve the
        % plot with a cosine?
        offsetcors=acmat(:,centercol);
        
        % find max correlation
        [maxcorr,offsetind]=max(offsetcors);
        % and the offset in time
        offsettime(i)=actimebins(offsetind);
        
        % now change our newts by that offset
        newts(i)=newts(i)-offsettime(i);
        corvals(i)=maxcorr;
        
        % get the NEW timebins for this event
        trialbints1=bins+newts(i);
        % revise our mean curve matrix
        allcurves(:,:,i)=makecurve(trialbints1,unitdata,binsz);
        
        waitbar(i/length(ts))
        
    end
close(h)
end

%
% figure; subplot(2,4,1); imagesc(smooth1s,'Ydata',bins);
% xlabel('cell number'); ylabel('time in seconds');
% subplot(2,4,2); imagesc(smooth2s,'Ydata',bins);
% xlabel('cell number'); ylabel('time in seconds');
% subplot(2,4,3:4);
% imagesc(acmat);
% ylabel('time in ms increments');
% subplot(2,4,5:8);
% plot(actimebins,acmat(:,centercol));
% hold on;
% plot(actimebins(offsetind),maxcorr,'r*');
%
%
% title(['max corr at ' num2str(offsettime) ' milliseconds']);
% xlabel(' <- neg offset in sec   pos offset in sec ->');

end
%% the curve finding function

function [smoothed]=makecurve(timestamps,unitdata,binsize)

% get the tuning curve for the event

[curves2]=eventspikematrix(timestamps(:),unitdata,0,binsize);
% get back orginal spike number
curves2=curves2.*binsize;
% kill cells that spike fewer times than 20
rates=sum(curves2);
%curves2(:,rates<10)=nan;
% smooth each by 5 pixels, or about 25 mSec
for i=1:size(curves2,2)
    smooth2(:,i)=smooth(curves2(:,i),9);
end
% zscore so no curve dominates cor matrix
smoothed=zscore(smooth2);
%smoothed=smooth2;
end
