function graphTimeCorrData(timecor,timeseries)

% graph the dimension x time corr keene


% retrofitted for imbedded struct data
CPI=timecor.CPI;
CPxI=timecor.CPxI;
CP=timecor.CP;
CxP=timecor.CxP;
xC=timecor.xC;
CxPI=timecor.CxPI;
CxPxI=timecor.CxPxI;

stats=0;
%% This is context coding for your region
figure;
subplot(2,2,1);
set(gca,'fontsize',16)

y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),CxP,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),CxP,'uni',0)),1);

xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[1 .8 0],'EdgeColor','none')
pl(1) = line(x,y,'linewidth',2,'color','r');

hold on


y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),xC,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),xC,'uni',0)),1);

xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[0 .8 1],'EdgeColor','none')
pl(2) = line(x,y,'linewidth',2,'color','b');


xlabel('Time to sample (s)')
ylabel('Correlation coefficient (r)')
axis tight
legend(pl,'Same Context','Diff. Context')

%is the observed t-stat at every point higher than the threshold (pval) for the
%bootstrap? with a very stringent p value

pval = .001;

over=cell2mat(CxP);  under=cell2mat(xC);
for h=1:length(under(1,:))
    cxt(h)=dprime(over(:,h),under(:,h));
end
% also taking out a,b,c and making them dummys
if stats==1

    h1 = double(cxt>=prctile(boot.ses.Cd,99.9));
    h1(h1==0) = nan;
    
    plot(timeseries,.23*h1,'k','linewidth',5)
end
ylim([-.1 .3])



%% This is position coding for your region
subplot(2,2,2);
set(gca,'fontsize',16)

y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),CP,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),CP,'uni',0)),1);

xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[1 .8 0],'EdgeColor','none')
pl(1) = line(x,y,'linewidth',2,'color','r');

hold on


y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),CxP,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),CxP,'uni',0)),1);

xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[0 .8 1],'EdgeColor','none')
pl(2) = line(x,y,'linewidth',2,'color','b');


xlabel('Time to sample (s)')
ylabel('Correlation coefficient (r)')
axis tight
legend(pl,'Same position','Diff. position')

%is the observed t-stat at every point higher than the threshold (pval) for the
%bootstrap? with a very stringent p value

pval = .001;


over=cell2mat(CP);  under=cell2mat(CxP);
for h=1:length(under(1,:))
    pos(h)=dprime(over(:,h),under(:,h));
end
% also taking out a,b,c and making them dummys
if stats==1

    h1 = double(pos>=prctile(boot.ses.Pd,99.9));
    h1(h1==0) = nan;
    
    plot(timeseries,.23*h1,'k','linewidth',5)
end
ylim([-.1 .3])



%% this is item in position coding for your region


subplot(2,2,4);
set(gca,'fontsize',16)

y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),CPI,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),CPI,'uni',0)),1);

% get shape of patch
xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[1 .8 0],'EdgeColor','none')
pl(1) = line(x,y,'linewidth',2,'color','r');



hold on


y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),CPxI,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),CPxI,'uni',0)),1);

xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[0 .8 1],'EdgeColor','none')
pl(2) = line(x,y,'linewidth',2,'color','b');


xlabel('Time to sample (s)')
ylabel('Correlation coefficient (r)')
axis tight
legend(pl,'Same item in place','Diff. item in place')
%jay edit those used to be a,b,c
% this is to get the significance, ttest because its bootstrapped
% i think this should be paired

over=cell2mat(CPI);  under=cell2mat(CPxI);
for h=1:length(under(1,:))
    itpos(h)=dprime(over(:,h),under(:,h));
end
if stats==1
    h = double(itpos>=prctile(boot.ses.Id,99.9)); % top .1%
    h(h==0) = nan;
    
    plot(timeseries,.23*h,'k','linewidth',5)
    
    ylim([-.1 .3])
end

%% this is item coding for your region


subplot(2,2,3);
set(gca,'fontsize',16)

y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),CxPI,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),CxPI,'uni',0)),1);

% get shape of patch
xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[1 .8 0],'EdgeColor','none')
pl(1) = line(x,y,'linewidth',2,'color','r');



hold on


y=nanmean(cell2mat(cellfun(@(a) nanmean(a,1),CxPxI,'uni',0)),1);
x=timeseries;
z=SEM(cell2mat(cellfun(@(a) nanmean(a,1),CxPxI,'uni',0)),1);

xx=[x x(end:-1:1)];
yy=[y-z y(end:-1:1)+z(end:-1:1)];
p=patch(xx,yy,'w');

set(p,'FaceColor',[0 .8 1],'EdgeColor','none')
pl(2) = line(x,y,'linewidth',2,'color','b');


xlabel('Time to sample (s)')
ylabel('Correlation coefficient (r)')
axis tight
legend(pl,'Same item','Diff. item')
%jay edit those used to be a,b,c
% this is to get the significance, ttest because its bootstrapped
% i think this should be paired


over=cell2mat(CxPI);  under=cell2mat(CxPxI);
for h=1:length(under(1,:))
    item(h)=dprime(over(:,h),under(:,h));
end
if stats==1
    h = double(item>=prctile(boot.ses.Id,99.9)); % top .1%
    h(h==0) = nan;
    
    plot(timeseries,.23*h,'k','linewidth',5)
    
    ylim([-.1 .3])
end


%% Now Item vs position coding, the lines above are where its >0

% timeseries is x axis
% pos is the y axis for pos

figure;



hold on
%
%plot(timeseries,cxt,'g','linewidth',3)
plot(timeseries,pos,'r','linewidth',3)
plot(timeseries,item,'b','linewidth',3)
plot(timeseries,itpos,'c','linewidth',3)

if stats==1
plot(timeseries,h*52,'linewidth',10,'color','b')
plot(timeseries,h1*55,'linewidth',10,'color','r')
end
set(gca,'fontsize',16)
ylabel('dprime')
xlabel('Time to sample (s)')
xlim([-3 3])
legend('Position','Item','itemXpos')
axis tight

end