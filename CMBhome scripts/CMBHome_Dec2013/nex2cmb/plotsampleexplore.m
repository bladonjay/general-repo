function [exploregrid,samplegrid] = plotsampleexplore(data,scale, gausswin, gaussdev)
%plotsampleexplore plots two heatplots under a common scale
% the first heatplot is for the object exploration period, bracketed by the
% trial start and board down flags.  The second heatplot is the object
% sample period bracketed by divider up and the choicepoints


% could add input parser but thats a pain... but input args is getting long
if ~exist('scale','var')
    data.spatial_scale=.5;
elseif isempty(scale)
    data.spatial_scale=.5;
else
    data.spatial_scale=scale;
end

if ~exist('gausswin','var')
    gausswin=5;
elseif isempty(gausswin)
    gausswin=5;
end

if ~exist('gaussdev','var')
    gaussdev=1;
elseif isempty(gaussdev)
    gaussdev=1;
end




% probably should add u input argument to determine how many, or which
% cells to plot.

for i=1:length(data.spike)
    figure(i);
    startevent=data.Label2Time('Trial Start');
    stopevent=data.Label2Time('Divider Down');
    data.epoch= [startevent stopevent];
    
    [expratemap]=data.plot_rate_map2([1,i],1,[],gausswin,gaussdev);
    ymax=max(max(expratemap));
    
    startevent=data.Label2Time('Divider Up');
    stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
    %stopevent=startevent+2;
    %startevent=stopevent-2;
    data.epoch= [startevent stopevent];
    [sampratemap]=plot_rate_map2(data,[1,i],1,[],gausswin,gaussdev);
    ymax=[ymax max(max(sampratemap))];
    % plot the first graph;
    subplot(2,1,1);
    startevent=data.Label2Time('Trial Start');
    stopevent=data.Label2Time('Divider Down');
    data.epoch= [startevent stopevent];
    plot_rate_map2(data,[1,i],0,[0 max(ymax)],gausswin,gaussdev);
    % plot second graph;
    subplot(2,1,2);
    startevent=data.Label2Time('Divider Up');
    stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
    %stopevent=startevent+2;
    %startevent=stopevent-2;
    data.epoch= [startevent stopevent];
    data.plot_rate_map2([1,i],0,[0 max(ymax)],gausswin,gaussdev);
end

end

