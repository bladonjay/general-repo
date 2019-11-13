function  plotallspatials(data,Session)
%This should plot a series of 6 spatial plots, 1 is exploration 2 is x on
%top, 3 is open field, 4 is all samples, 5 is y on top and 6 should be the
%spatial ac.  I'm hoping they all go on the same scale(at least the four
%maze plots)
for i=1:length(data.spike)
    figure(i);
    
    startevent=data.Label2Time('Trial Start');
    stopevent=data.Label2Time('Divider Down');
    data.epoch= [startevent stopevent];
    
    [expratemap]=data.plot_rate_map([1,i],1);
    ymax=max(max(expratemap));
    
    startevent=data.Label2Time('Divider Up');
    stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
   
    
    data.epoch= [startevent stopevent];
    [sampratemap]=data.plot_rate_map([1,i],1);
    ymax=[ymax max(max(sampratemap))];
    
    % plot the first graph;
    sp(1)=subplot(2,3,1);
    startevent=data.Label2Time('Trial Start');
    stopevent=data.Label2Time('Divider Down');
    data.epoch= [startevent stopevent];
    data.plot_rate_map([1,i],0,[],[],[0 max(ymax)]);
    
    % plot the second graph;
    % USING 1 second before and 3 after, although I probably could code in
    % board to eat if i really wanted to
    sp(2)=subplot(2,3,2);
    topsamples=sort([Session.bf.XLE Session.bf.YRE Session.bf.XRW Session.bf.YLW]);
    topepoch=transpose([topsamples-1;topsamples+3]);
    data.epoch=topepoch;
    data.plot_rate_map([1,i]);
    
    
    % Plot third graph;
    sp(3)=subplot(2,3,3);
    tempstartevent=data.Label2Time('Open Field Flags');
    data.epoch=[tempstartevent(1) tempstartevent(end)];
    data.plot_rate_map([1,i]);
   
    
    % plot fourth graph;
    sp(4)=subplot(2,3,4);
    startevent=data.Label2Time('Divider Up');
    stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
    data.epoch= [startevent stopevent];
    data.plot_rate_map([1,i]);
    
    % plot fifth graph;
    sp(5)=subplot(2,3,5);
    bottomsamples=sort([Session.bf.XRE Session.bf.YLE Session.bf.XLW Session.bf.YRW]);
    bottomepoch=transpose([bottomsamples-1;bottomsamples+3]);
    data.epoch=bottomepoch;
    data.plot_rate_map([1,i]);
    
    % plot sixth graph
    sp(6)=subplot(2,3,6);
    tempstartevent=data.Label2Time('Open Field Flags');
    data.epoch=[tempstartevent(1) tempstartevent(end)];
    data.plot_rate_map_ac([1,i]);
    
    
end