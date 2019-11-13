function [fighandle,filelocation] = RainbowRaster(samples,unitdata,unitnumber,saveout)
% Make a pretty raster plot




fighandle=0; filelocation=[];
if nargin<4
    saveout=0;
end
binnums=20;
colors=parula(7);
colors=colors(1:6,:);
un=unitnumber;
secbefore=2; secafter=3;


figure('Position',[-999   217   850   687]);
try
    
    % gotta get left, right, beads and straws
    beadsleft =  sum(samples(:,[2 3])==[1 1],2)==2;
    beadsright = sum(samples(:,[2 3])==[2 1],2)==2;
    strawsleft = sum(samples(:,[2 3])==[1 2],2)==2;
    strawsright= sum(samples(:,[2 3])==[2 2],2)==2;
    
    subplot(2,8,1:2); thispos=get(gca,'Position');
    [a,b,c,d,e]=EventPethraPlot(unitdata.units(un).ts,samples(beadsleft,1),'AxisCoords',thispos,...
        'EventID', samples(beadsleft,6),'EventSort',[1:6],'SecBefore',secbefore,...
        'SecAfter',secafter,'EventColors',colors,'histbins',binnums);
    [a,b,c,d]=event_spikes(unitdata.units(un).ts,samples(beadsleft,1),.5,.5);
    plots=get(gcf,'children');
    topplot=get(plots(2),'Position'); oldytick=get(plots(2),'YTick');
    axes('Position',topplot+[.175 0 -.08 0]);
    %b=smooth(b);
    for it=1:length(b)-1
        plot([b(it) b(it+1)],[it+.5 it+1.5],'Color',e(it,5:7),'LineWidth',2);
        hold on
    end
    ylim([1 length(b)*1]);
    
    
    subplot(2,8,5:6); thispos=get(gca,'Position');
    [a,b,c,d,e]=EventPethraPlot(unitdata.units(un).ts,samples(beadsright,1),'AxisCoords',thispos,...
        'EventID', samples(beadsright,6),'EventSort',[1:6],'SecBefore',secbefore,...
        'SecAfter',secafter,'EventColors',colors,'histbins',binnums);
    [a,b,c,d]=event_spikes(unitdata.units(un).ts,samples(beadsright,1),.5,.5);
    plots=get(gcf,'children');
    topplot=get(plots(2),'Position'); oldytick=get(plots(2),'YTick');
    axes('Position',topplot+[.175 0 -.08 0]);
    %b=smooth(b);
    for it=1:length(b)-1
        plot([b(it) b(it+1)],[it+.5 it+1.5],'Color',e(it,5:7),'LineWidth',2);
        hold on
    end
    ylim([1 length(b)+1]);
    
    
    subplot(2,8,9:10); thispos=get(gca,'Position');
    [a,b,c,d,e]=EventPethraPlot(unitdata.units(un).ts,samples(strawsleft,1),'AxisCoords',thispos,...
        'EventID', samples(strawsleft,6),'EventSort',[1:6],'SecBefore',secbefore,...
        'SecAfter',secafter,'EventColors',colors,'histbins',binnums);
    [a,b,c,d]=event_spikes(unitdata.units(un).ts,samples(strawsleft,1),.5,.5);
    plots=get(gcf,'children');
    topplot=get(plots(2),'Position'); oldytick=get(plots(2),'YTick');
    axes('Position',topplot+[.175 0 -.08 0]);
    %b=smooth(b);
    for it=1:length(b)-1
        plot([b(it) b(it+1)],[it+.5 it+1.5],'Color',e(it,5:7),'LineWidth',2);
        hold on
    end
    ylim([1 length(b)+1]);
    
    
    subplot(2,8,13:14); thispos=get(gca,'Position');
    [a,b,c,d,e]=EventPethraPlot(unitdata.units(un).ts,samples(strawsright,1),'AxisCoords',thispos,...
        'EventID', samples(strawsright,6),'EventSort',[1:6],'SecBefore',secbefore,...
        'SecAfter',secafter,'EventColors',colors,'histbins',binnums);
    [a,b,c,d]=event_spikes(unitdata.units(un).ts,samples(strawsright,1),.5,.5);
    plots=get(gcf,'children');
    topplot=get(plots(2),'Position'); oldytick=get(plots(2),'YTick');
    axes('Position',topplot+[.175 0 -.08 0]);
    %b=smooth(b);
    for it=1:length(b)-1
        plot([b(it) b(it+1)],[it+.5 it+1.5],'Color',e(it,5:7),'LineWidth',2);
        hold on
    end
    ylim([1 length(b)+1]);
    
    
    %thisname=[FinishedRat(ses).name ' ' FinishedRat(ses).unitdata.units(un).name];
    %thisname(thisname=='_')='-';
    %titlebox= uicontrol(gcf,'Style','Text','String',thisname,...
    %    'Position',[300 640 200 30],'units','pixels');
    
    %text(gcf,.5,.5,thisname,'Units','normalized');
    smallplots=get(gcf,'Children');
    linkaxes(smallplots(2:3:length(smallplots)));
    linkaxes(smallplots(1:3:length(smallplots)),'x');
    % this will link all x axes of the side plots, and then stretch them
    xlimits=get(smallplots(1),'Xlim'); xlimits(2)=xlimits(2)*1.15;
    

    sptypes=[repmat(1:3,1,4)];
    for sps=1:length(smallplots)
        set(smallplots(sps),'FontSize',12);
        if sptypes(sps)==1
            % righthand plot
            set(smallplots(sps).XLabel,'String','Rate');
            set(smallplots(sps),'YTickLabel',[],'XDir',...
                'Reverse','YDir','Reverse');
            
        elseif sptypes(sps)==2
            % bottom plot
            set(smallplots(sps).XLabel,'String','Time');
            set(smallplots(sps).YLabel,'String','Rate');
        elseif sptypes(sps)==3
            set(smallplots(sps).YLabel,'String','Trial Number');
        end
    end
    
    set(smallplots(1),'xlim',xlimits);
    set(smallplots(2),'Xlim',[-secbefore secafter])
    set(smallplots(3),'Xlim',[-secbefore secafter])
    fighandle=gcf;
    if saveout
        savefig(gcf,[savedir dash thisname]);
        saveas(gcf,[savedir dash thisname],'jpg');
        close
        filelocation=savedir;
    end
catch
    fprintf('something got fucked up \n');
    close
    filelocation=0;
    fighandle=0;
end
end

