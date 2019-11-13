% THis is a SOP for generating a variety of spatial plots using CMBHOME
% code
%%%%%%%

%have folder on hand to analyze sessions (folder of nex files)

[data, nex] = jayNEX2CMBbatch('',0);

% This will first ask you first for the directory of nex files, then ask
% you to choose the mat file with the corrected xy coordinate data

%%

uiopen
%%

% This figure generator creates paired figures- top is context exploration
% period(you can modify it to be just the first ten seconds if youd like)
% the bottom plot is the item sampling epoch from board up to start eat

plotsampleexplore(data,.5,10,1);


%% To plot all trials from trial start to decision point

for i=1:length(data.spike)
    figure(i*2);
    
    startevent=data.Label2Time('Trial Start');
    stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
    data.epoch= [startevent stopevent];
    data.plot_rate_map2([1,i]);
    
    
end


%% Lets see if we can get one context only, and plot the ac for that side:


% Lets start with the right context
startevent=data.Label2Time('Right Context Enter');
stopevent=data.Label2Time('Right Context Exit');
% have to cut it off at 90 because sometimes entrys and exits arent equal
% (last trial)

shorter=min([length(startevent) length(stopevent)]);
startevent=startevent(1:shorter);
stopevent=stopevent(1:shorter);
data.epoch=[startevent stopevent];
for i=3:length(data.spike)
    figure(i);
    subplot(1,2,2);
    data.plot_rate_map([1,i]);
    %[gridness, props]  = data.Gridness([1,i]);
    %cellprops(i).right=props(1);
    title('right context only');
end
% now we plot left context
startevent=data.Label2Time('Left Context Enter');
stopevent=data.Label2Time('Left Context Exit');
% have to cut it off at 90 because sometimes entrys and exits arent equal
% (last trial)
shorter=min([length(startevent) length(stopevent)]);
startevent=startevent(1:shorter);
stopevent=stopevent(1:shorter);
data.epoch=[startevent stopevent];

data.epoch=[startevent stopevent];
for i=3:length(data.spike)
    figure(i);
    subplot(1,2,1);
    data.plot_rate_map([1,i]);
    %[gridness, props]  = data.Gridness([1,i]);
    %cellprops(i).left=props(1);
    title('left context only');
end

%% set to open field
% there clearly needs to be a velocity filter on this...
ofdata=data.SetEpoch;
for i=1:length(data.spike)
    figure(i);
    subplot(1,2,1);
    ofdata.plot_rate_map2([1,i]);
    subplot(1,2,2);
    ofdata.plot_rate_map_ac([1,i]);
end

%% Gridness scores and grid metrics;

[gridness,props]=ofdata.Gridness([1,1]);

% Also
% spatial_information = ofdata.SpatialInformation([1,1]);

%%  To look at first 30 vs last 30 trials
for i=1:length(data.spike)
    figure(i);
cellnumber=[1,i];

startevent=data.Label2Time('Divider Up');
stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
startevent=startevent(1:30); stopevent=stopevent(1:30);
data.epoch= [startevent stopevent];
subplot(2,1,1);
data.plot_rate_map(cellnumber);

startevent=data.Label2Time('Divider Up');
stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
startevent=startevent(60:end); stopevent=stopevent(60:end);
data.epoch= [startevent stopevent];
subplot(2,1,2);
data.plot_rate_map(cellnumber);
end
%%
% This is to create the two maps in komo (top xy in one config, bottom xy
% in a different config)
% I hardcoded in the sample times, but sams matrix would be better for
% sample times- you can use his sample lengths from 'samples' in his saved
% .mats

%  YOU MUST CHOOSE THE CORRECT SESSION MAT FILE
uiopen
topsamples=sort([Session.bf.XLE Session.bf.YRE Session.bf.XRW Session.bf.YLW]);
bottomsamples=sort([Session.bf.XRE Session.bf.YLE Session.bf.XLW Session.bf.YRW]);

topepoch=transpose([topsamples-1;topsamples+2]);
bottomepoch=transpose([bottomsamples-1;bottomsamples+2]);
for i=1:length(data.spike)
    figure(i);
    subplot(2,1,1);
    data.epoch=topepoch;
    data.plot_rate_map([1,i]);
    
    subplot(2,1,2);
    data.epoch=bottomepoch;
    data.plot_rate_map([1,i]);
end

%%  Maybe I can use the CMB epoching to create some nonspatial scripts and commands
% waste of time...

% First I'll see if i can get a generalized pSPIKE | event time
data=data.SetEpoch;
% get spikes and time for whole session

tottime=data.ts(end)-data.ts(1);
for i=1:length(data.spike)
    totspikes(i)=length(data.spike(i).ts);
end

% this is for the object sampling part
startevent=(data.Label2Time('Divider Up'));
stopevent=sort([data.Label2Time('Correct Choice');data.Label2Time('Incorrect Choice')]);
%stopevent=startevent+2;
%startevent=stopevent-2;
data.epoch= [startevent stopevent];
% now get spikes and time for epoch

epochtime=sum(diff(data.epoch,1,2));
for i=1:length(data.spike)
epochspikes(i)=length(data.spk_ts([1,i]));
end

% now compare ratios
p_epoch=epochtime/tottime;
p_epspikes=epochspikes/totspikes;


spikeratio=p_epspikes/p_epoch;  % are spikes in the epoch more or less likely than overall?

% 1=equally likely, more than 1 more likely, less than 1 less likely


%data.spk_ts([1,1])

