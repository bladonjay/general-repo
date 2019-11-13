function [ax]=PrepAxes(inputax)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(inputax)
    if contains(get(inputax(i),'type'),'axes')
        % gather all its information
        myplot=getframe(inputax(i));
        myimage=frame2im(myplot);
        mypos=get(inputax(i),'Position');
        % now axes
        if ~isempty(get(inputax(i),'XTickLabel'))
            oldlimx=get(inputax(i),'XLim');
            oldtickx=get(inputax(i),'XTick');
            oldlabelx=get(inputax(i),'XTickLabel');
            newlimx=size(myimage,2)-1;
            stolenlabelx=get(inputax(i).XLabel,'String');
            if strcmpi(get(inputax(i),'XDir'),'reverse')
                newtickx=fliplr(newlimx-(oldtickx-oldlimx(1))...
                    .*(newlimx/(oldlimx(2)-oldlimx(1))))+1;
                oldlabelx=flipud(oldlabelx);
            else
                newtickx=(oldtickx-oldlimx(1)).*(newlimx/(oldlimx(2)-oldlimx(1)))+1;
            end
        end
        if ~isempty(get(inputax(i),'YTickLabel'))
            oldlimy=get(inputax(i),'YLim');
            oldticky=get(inputax(i),'YTick');
            oldlabely=get(inputax(i),'YTickLabel');
            newlimy=size(myimage,1)-1;
            stolenlabely=get(inputax(i).YLabel,'String');
            stolentitle=get(inputax(i),'Title');
            if ~strcmpi(get(inputax(i),'YDir'),'reverse')
                newticky=fliplr(newlimy-(oldticky-oldlimy(1))...
                    .*(newlimy/(oldlimy(2)-oldlimy(1))))+1;
                oldlabely=flipud(oldlabely);
            else
                newticky=(oldticky-oldlimy(1)).*(newlimy/(oldlimy(2)-oldlimy(1)))+1;
            end
        end
        
        % now make it so that we replace the old axis
        
        % now replace the current axes and clear its ticks
        cla(inputax(i),'reset'); set(inputax(i),'XTick',[],'YTick',[]);
        ha=axes('Position',mypos); 
        imagesc(myimage);
        
        % add x labels
        if exist('oldlimx','var')
            set(ha,'XTick',newtickx,'XTickLabel',oldlabelx);
            try set(ha.XLabel,'String',stolenlabelx); end
        else
            set(ha,'XTick',[]);
        end
        
        
        % add ylabels
        if exist('oldlimy','var')
            set(ha,'YTick',newticky,'YTickLabel',oldlabely);
            try set(ha.YLabel,'string',stolenlabely); end
        else
            set(ha,'YTick',[]);
        end
        
    end
end


end

