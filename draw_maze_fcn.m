function [xreturn,yreturn,track] = draw_maze_fcn(guihandle,fighandle,nexthandle,donehandle)
	set(fighandle,'windowbuttondownfcn',@wbdcb);
	%set up the drawing mode
	xvals = [];
	yvals = [];
	xreturn = [];
	yreturn = [];
    track = [];
    color='mbyrcgk';
    %
    %also start a new series of lines for each new track.
    
	function wbdcb(src,evnt) %Adds clicked point to list of x-,y-values
     if strcmp(get(src,'SelectionType'),'normal')
        set(src,'pointer','fullcrosshair')
        cp = get(gca,'CurrentPoint');
        xinit = cp(1,1);
        yinit = cp(1,2);
        m=rem(get(nexthandle, 'UserData'),7)+1;
        xvals = [xvals, xinit];
        yvals = [yvals, yinit];
        track = [track, get(nexthandle, 'UserData')];
        h1 = line('XData',xinit,'YData',yinit,'Marker','o','color',color(m),'LineWidth',3);
        set(src,'WindowButtonMotionFcn',@wbmcb)
        set(src,'WindowButtonUpFcn',@wbucb)
     end

        function wbmcb(src,evnt) %Temporarily plots line while user moves 
            %towards the next point
           cp = get(gca,'CurrentPoint');
           xdat = [xinit,cp(1,1)];
           ydat = [yinit,cp(1,2)];
           m=rem(get(nexthandle, 'UserData'),7)+1;
           set(h1,'XData',xdat,'YData',ydat,'color',color(m));
           drawnow
        end

        function wbucb(src,evnt)
           if strcmp(get(src,'SelectionType'),'alt')
              cp = get(gca,'CurrentPoint');
			  xvals = [xvals,cp(1,1)];
			  yvals = [yvals,cp(1,2)];
              curtrack=get(nexthandle, 'UserData');
              track = [track, curtrack];
              set(nexthandle, 'UserData',curtrack+1);%Advance to next track
              set(src,'Pointer','arrow');
              set(src,'WindowButtonMotionFcn','');
              set(src,'WindowButtonUpFcn','');
			  xreturn = xvals';
			  yreturn = yvals';
           else
               return
           end
        end
    end
waitfor(donehandle,'Value',1);
close(guihandle);

xreturn = xvals';
yreturn = yvals';
end
