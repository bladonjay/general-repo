function [fig]=learningrasterplot(session, unit)


enters1=session.bf.r_enter'; enters1(:,2)=1;
enters2=session.bf.l_enter'; enters2(:,2)=2;
enters=sortrows([enters1; enters2],1); enters(:,3)=1:length(enters);

clear trialthirds
trialthirds(:,1)=session.samples(:,11)<45;
trialthirds(:,2)=session.samples(:,11)>45;

samples=session.samples;

% this is sloppy gonna have to fix it
j=unit;
fig=figure;

% probably can clean this up too
rowpos=[1 2 3 4 5 6 7; 8 9 10 11 12 13 14];
for k=1:2
    % left enters
    ax=subplot(2,7,rowpos(k,1)); axcoords=get(ax,'Position');
    leftenters=enters(enters(:,2)==2 & ceil((enters(:,3)/size(enters,1))*3)==k,:);
    EventPethraPlot(session.units(j).ts',leftenters(:,1),...
        'AxisCoords',axcoords);
    xlabel('Left Enters','FontSize',12);
    % right enters
    ax=subplot(2,7,rowpos(k,2)); axcoords=get(ax,'Position');
    rightenters=enters(enters(:,2)==1 & ceil((enters(:,3)/size(enters,1))*3)==k,:);
    EventPethraPlot(session.units(j).ts',rightenters(:,1),...
        'AxisCoords',axcoords);
    xlabel('Right Enters','FontSize',12);
    
    % top left third col is 1, 5th col is 1?
    ax=subplot(2,7,rowpos(k,4)); axcoords=get(ax,'Position');
    thesetrials=samples(trialthirds(:,k)==1 & samples(:,3)==1 & samples(:,5)==1,:);
    EventPethraPlot(session.units(j).ts',thesetrials,'EventID',thesetrials(:,7),...
        'AxisCoords',axcoords,'EventSort',[1 2],'EventColor',[0 0 1; 1 0 0]);
    xlabel('Left, West','FontSize',12);
    
    % bottom left third col is 2, 5th col is 1?
    ax=subplot(2,7,rowpos(k,5)); axcoords=get(ax,'Position');
    thesetrials=samples(trialthirds(:,k)==1 & samples(:,3)==2 & samples(:,5)==1,:);
    EventPethraPlot(session.units(j).ts',thesetrials,'EventID',thesetrials(:,7),...
        'AxisCoords',axcoords,'EventSort',[1 2],'EventColor',[0 0 1; 1 0 0]);
    xlabel('Right, West','FontSize',12);
    
    % top right third col is 1, 5th col is 2?
    ax=subplot(2,7,rowpos(k,6)); axcoords=get(ax,'Position');
    thesetrials=samples(trialthirds(:,k)==1 & samples(:,3)==2 & samples(:,5)==1,:);
    EventPethraPlot(session.units(j).ts',thesetrials,'EventID',thesetrials(:,7),...
        'AxisCoords',axcoords,'EventSort',[1 2],'EventColor',[0 0 1; 1 0 0]);
    xlabel('Left East','FontSize',12);
    
    % top right third col is 2, 5th col is 2?
    ax=subplot(2,7,rowpos(k,7)); axcoords=get(ax,'Position');
    thesetrials=samples(trialthirds(:,k)==1 & samples(:,3)==2 & samples(:,5)==2,:);
    EventPethraPlot(session.units(j).ts',thesetrials,'EventID',thesetrials(:,7),...
        'AxisCoords',axcoords,'EventSort',[1 2],'EventColor',[0 0 1; 1 0 0]);
   xlabel('Right East','FontSize',12);

    
end

allaxes=get(gcf,'Children');
linkaxes(allaxes(1:2:end));


end


