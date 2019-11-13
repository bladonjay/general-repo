function [coords, events ] = PlotFlagCoords(LEDdata, events, eventtypes)
% function [coords, events ] = PlotFlagCoords(unitdata, samples)
% inputs: 
%   LEDdata: must be Mx3 or Mx5 matrix, first col is always ts, second
%   and third are led1, and 4 and 5 are led2.
%   events: must be column of timestamps for indexing within timestamps
%   eventtypes: idetifier for the classes of events, ought to be four.

after=.1; before=-.1;

figure
% find the ts within
X1=cellfun(@(a) nanmean(LEDdata((a+after>LEDdata(:,1) & a+before< LEDdata(:,1)),2)),num2cell(events(:,1),2));
Y1=cellfun(@(a) nanmean(LEDdata((a+after>LEDdata(:,1) & a+before< LEDdata(:,1)),3)),num2cell(events(:,1),2));
%if size
X2=cellfun(@(a) nanmean(LEDdata((a+after>LEDdata(:,1) & a+before< LEDdata(:,1)),4)),num2cell(events(:,1),2));
Y2=cellfun(@(a) nanmean(LEDdata((a+after>LEDdata(:,1) & a+before< LEDdata(:,1)),5)),num2cell(events(:,1),2));
X=nanmean([X1 X2],2);
Y=nanmean([Y1 Y2],2);
coords=[X Y];


% there are four poitions, plot each of them in the positions according to
% the tracking
plot(coords(eventtypes==4,1),coords(eventtypes==4,2),'go')
hold on
plot(coords(eventtypes==1,1),coords(eventtypes==1,2),'x')
plot(coords(eventtypes==2,1),coords(eventtypes==2,2),'rx')
plot(coords(eventtypes==3,1),coords(eventtypes==3,2),'ko')

% and this tracks the timestamps so when you click on a datapoint you can
% grab its ts
CheckPos(X,Y,events(:,1))


end

