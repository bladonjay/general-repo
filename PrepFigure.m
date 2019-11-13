
function [outfig]=PrepFigure(infig,figtype)
% PrepFigure

% function PrepFigure
% this just takes all the plots and turns them into pictures, it preserves
% the axes so they can be imported as text


%%
if nargin<1
    hfig=get(gcf);
else
    hfig=infig;
end
mychildren=get(gcf,'Children');

if ~exist('figtype','var')
    figtype=1;
elseif ~isnumeric(figtype)
    figtype=1;
end

% normal figure
if figtype==1
figure('position',hfig.Position);
for i=1:length(mychildren)
    if contains(mychildren(i).Type,'axes')
        myplot=getframe(mychildren(i));
        myimage=frame2im(myplot);
        axes('Position',get(mychildren(i),'Position'));
        % here i'm leaving the axis reversed, probably is more complicated, but
        % whatever. an alternative approach is to plot my image with the
        % scaling baked in and set the ticks that way
        imagesc(myimage);
        set(gca,'FontSize',get(mychildren(i),'FontSize'));
        if ~isempty(get(mychildren(i),'XTickLabel'))
            oldlim=get(mychildren(i),'XLim');
            oldtick=get(mychildren(i),'XTick');
            oldlabel=get(mychildren(i),'XTickLabel');
            newlim=size(myimage,2)-1;
            if strcmpi(get(mychildren(i),'XDir'),'reverse')
                newtick=fliplr(newlim-(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1))))+1;
                set(gca,'XTick',newtick,'XTickLabel',flipud(oldlabel));
            else
                newtick=(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1)))+1;
                set(gca,'XTick',newtick,'XTickLabel',oldlabel);
            end
        else
            set(gca,'XTick',[]);
        end
        
        if ~isempty(get(mychildren(i),'YTickLabel'))
            oldlim=get(mychildren(i),'YLim');
            oldtick=get(mychildren(i),'YTick');
            oldlabel=get(mychildren(i),'YTickLabel');
            newlim=size(myimage,1)-1;
            if ~strcmpi(get(mychildren(i),'YDir'),'reverse')
                newtick=fliplr(newlim-(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1))))+1;
                set(gca,'YTick',newtick,'YTickLabel',flipud(oldlabel));
            else
                newtick=(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1)))+1;
                set(gca,'YTick',newtick,'YTickLabel',oldlabel);
            end
        else
            set(gca,'YTick',[]);
        end
        stolenlabel=get(mychildren(i),'XLabel'); set(gca,'XLabel',stolenlabel);
        stolenlabel=get(mychildren(i),'YLabel'); set(gca,'YLabel',stolenlabel);
        set(gca,'Title',get(mychildren(i),'Title'));
    
    elseif contains(mychildren(i).Type,'colorbar')
        cb=colorbar;
        oldpos=get(mychildren(i),'Position');
        cb.Position=oldpos;
        oldlimits=get(mychildren(i),'Limits');
        cb.Limits=oldlimits;
        oldticks=get(mychildren(i),'Ticks');
        cb.Ticks=oldticks;
        oldticklabels=get(mychildren(i),'TickLabels');
        cb.TickLabels=oldticklabels;
    end   
end
end

%% if 3d plots

if figtype==2
figure('position',hfig.Position);
for i=1:length(mychildren)
    if contains(mychildren(i).Type,'axes')
        myplot=getframe(mychildren(i));
        myimage=frame2im(myplot);
        axes('Position',get(mychildren(i),'Position'));
        % here i'm leaving the axis reversed, probably is more complicated, but
        % whatever. an alternative approach is to plot my image with the
        % scaling baked in and set the ticks that way
        imagesc(myimage);
        set(gca,'FontSize',get(mychildren(i),'FontSize'));
        if ~isempty(get(mychildren(i),'XTickLabel'))
            oldlim=get(mychildren(i),'XLim');
            oldtick=get(mychildren(i),'XTick');
            oldlabel=get(mychildren(i),'XTickLabel');
            newlim=size(myimage,2)-1;
            if strcmpi(get(mychildren(i),'XDir'),'reverse')
                newtick=fliplr(newlim-(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1))))+1;
                set(gca,'XTick',newtick,'XTickLabel',flipud(oldlabel));
            else
                newtick=(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1)))+1;
                set(gca,'XTick',newtick,'XTickLabel',oldlabel);
            end
        else
            set(gca,'XTick',[]);
        end
        % Y
        if ~isempty(get(mychildren(i),'YTickLabel'))
            oldlim=get(mychildren(i),'YLim');
            oldtick=get(mychildren(i),'YTick');
            oldlabel=get(mychildren(i),'YTickLabel');
            newlim=size(myimage,1)-1;
            if ~strcmpi(get(mychildren(i),'YDir'),'reverse')
                newtick=fliplr(newlim-(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1))))+1;
                set(gca,'YTick',newtick,'YTickLabel',flipud(oldlabel));
            else
                newtick=(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1)))+1;
                set(gca,'YTick',newtick,'YTickLabel',oldlabel);
            end
        else
            set(gca,'YTick',[]);
        end
        % Z is special...
        if ~isempty(get(mychildren(i),'ZTickLabel'))
            oldlim=get(mychildren(i),'ZLim');
            oldtick=get(mychildren(i),'ZTick');
            oldlabel=get(mychildren(i),'ZTickLabel');
            newlim=size(myimage,1)-1;
            if ~strcmpi(get(mychildren(i),'ZDir'),'reverse')
                newtick=fliplr(newlim-(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1))))+1;
                set(gca,'ZTick',newtick,'ZTickLabel',flipud(oldlabel));
            else
                newtick=(oldtick-oldlim(1)).*(newlim/(oldlim(2)-oldlim(1)))+1;
                set(gca,'ZTick',newtick,'ZTickLabel',oldlabel);
            end
        else
            set(gca,'ZTick',[]);
        end
        stolenlabel=get(mychildren(i),'ZLabel'); set(gca,'ZLabel',stolenlabel);
        stolenlabel=get(mychildren(i),'XLabel'); set(gca,'XLabel',stolenlabel);
        stolenlabel=get(mychildren(i),'YLabel'); set(gca,'YLabel',stolenlabel);
        set(gca,'Title',get(mychildren(i),'Title'));
    end
end
end


%% if other figtype
if figtype==5

set(gcf,'Position',[517   718   375   247]);
mychildren=get(gcf,'Children');
mypos=get(gcf,'Position');
figure('Position',mypos);
for i=1:2
        myplot=getframe(mychildren(i));
        oldlim=get(mychildren(i),'YLim');
        oldtick=get(mychildren(i),'YTick');
        oldticklabel=get(mychildren(i),'YTickLabel');
        myimage=frame2im(myplot);
        axes('Position',get(mychildren(i),'Position'));
        
        
        if i==2
        imagesc(myimage);
        % Convert the axes to a % of axis lim
        standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
        newlim=get(gca,'YLim');
        set(gca,'YTick',standardtick*newlim(2));
        set(gca,'YTickLabel',oldticklabel,'TickLength',[0 0]);
        
            %ylabel('Trial Number');
            set(gca,'XTick',[],'XTickLabel',[]);
            temp=get(mychildren(i),'Children');
            
            %ylabel('Trial Number');
        elseif i==1
            imagesc(flipud(myimage));
            standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
            newlim=get(gca,'YLim');
            set(gca,'YTick',standardtick*newlim(2));
            set(gca,'YTickLabel',oldticklabel,'TickLength',[0 0]);
            
            
            %ylabel('Rate (Hz)');
            set(gca,'ydir','normal');
            newlim=get(gca,'XLim');
            oldlim=get(mychildren(i),'XLim');
            oldtick=get(mychildren(i),'XTick');
            oldticklabel=get(mychildren(i),'XTickLabel');
            standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
            set(gca,'xTick',standardtick*newlim(2));
            set(gca,'xTickLabel',oldticklabel,'TickLength',[0 0]);
            temp=get(mychildren(i),'Children');
            
            %xlabel('Seconds From Delay Start'); ylabel('Rate (Hz)');
        end

end
end



%%
if figtype==1 && length(mychildren)>1
% contingency for 2 or 8 plots, this really just gets used for the event
% pethra plots anyways
mychildren=get(gcf,'Children');
mypos=get(gcf,'Position');
figure('Position',mypos);
if length(mychildren)==2
    for i=1:2
        myplot=getframe(mychildren(i));
        oldlim=get(mychildren(i),'YLim');
        oldtick=get(mychildren(i),'YTick');
        oldticklabel=get(mychildren(i),'YTickLabel');
        myimage=frame2im(myplot);
        axes('Position',get(mychildren(i),'Position'));
        
        
        if i==2
        imagesc(myimage);
        % Convert the axes to a % of axis lim
        standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
        newlim=get(gca,'YLim');
        set(gca,'YTick',standardtick*newlim(2));
        set(gca,'YTickLabel',oldticklabel,'TickLength',[0 0]);
        
            %ylabel('Trial Number');
            set(gca,'XTick',[],'XTickLabel',[]);
            temp=get(mychildren(i),'Children');
            try
                cellname=get(temp(1),'String');
                t=annotation('TextBox','Position',[0.1929 0.9310 0.6911 0.0500],...
                    'String',cellname,'LineStyle','none','HorizontalAlignment','center');
            catch
            end
            ylabel('Trial Number');
        elseif i==1
            imagesc(flipud(myimage));
            standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
            newlim=get(gca,'YLim');
            set(gca,'YTick',standardtick*newlim(2));
            set(gca,'YTickLabel',oldticklabel,'TickLength',[0 0]);
            
            
            ylabel('Rate (Hz)');
            set(gca,'ydir','normal');
            newlim=get(gca,'XLim');
            oldlim=get(mychildren(i),'XLim');
            oldtick=get(mychildren(i),'XTick');
            oldticklabel=get(mychildren(i),'XTickLabel');
            standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
            set(gca,'xTick',standardtick*newlim(2));
            set(gca,'xTickLabel',oldticklabel,'TickLength',[0 0]);
            temp=get(mychildren(i),'Children');
            try
                cellname=get(temp(1),'String');
                t=annotation('TextBox','Position',[0.1929 0.9310 0.6911 0.0500],...
                    'String',cellname,'LineStyle','none','HorizontalAlignment','center');
            catch
            end
            xlabel('Seconds From Delay Start'); ylabel('Rate (Hz)');
        end
    end
    
elseif length(mychildren)==8
    
    for i=1:length(mychildren)
        myplot=getframe(mychildren(i));
        oldlim=get(mychildren(i),'YLim');
        oldtick=get(mychildren(i),'YTick');
        oldticklabel=get(mychildren(i),'YTickLabel');
        myimage=frame2im(myplot);
        axes('Position',get(mychildren(i),'Position'));
        imagesc(flipud(myimage));
        set(gca,'ydir','normal');
        
        % Convert the axes to a % of axis lim
        standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
        newlim=get(gca,'YLim');
        set(gca,'YTick',standardtick*newlim(2));
        set(gca,'YTickLabel',oldticklabel,'TickLength',[0 0]);
        if mod(i,2)==0
            %ylabel('Trial Number');
            set(gca,'XTick',[],'XTickLabel',[]);
        else
            %ylabel('Rate (Hz)');
            newlim=get(gca,'XLim');
            oldlim=get(mychildren(i),'XLim');
            oldtick=get(mychildren(i),'XTick');
            oldticklabel=get(mychildren(i),'XTickLabel');
            standardtick=(oldtick-oldlim(1))/(oldlim(2)-oldlim(1));
            set(gca,'xTick',standardtick*newlim(2));
            set(gca,'xTickLabel',oldticklabel,'TickLength',[0 0]);
            %temp=get(mychildren(i),'Children');
            %cellname=get(temp(1),'String');
            %t=annotation('TextBox','Position',[0.1929 0.9310 0.6911 0.0500],...
            %'String',cellname,'LineStyle','none','HorizontalAlignment','center');
            xlabel('Seconds From Delay Start');
        end
    end
end
end