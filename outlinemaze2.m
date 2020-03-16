function [xpos,ypos,tracks] = outlinemaze2(xydata,im)
%OUTLINEMAZE Plots tracking data and allows user to outline the maze
%   xydata-2xN array of tracking data (x & y position)

%Check for proper format of xydata
if size(xydata,2)==3
    xydata(:,1)=[];
elseif size(xydata,2)~=2
    error('Tracking data must be a 2xN or 3xN matrix')
end

%Start GUI and grab handles
guihandle=trackgui;
guistructure=guidata(guihandle);
imagehandle=guistructure.axes1;
fighandle=guistructure.figure1;
nexthandle=guistructure.pushbutton1;
donehandle=guistructure.togglebutton2;

%Set image as backdrop
if exist('im','var')
    if ~isempty(im)
        imshow(im,'Parent',imagehandle, 'XData',[1 1024],'YData',[768 1]);
        hold(imagehandle,'on')
    end
end
%Plot tracking
plot(imagehandle,xydata(:,1),xydata(:,2),'Color',[0.5 0.5 0.5]);
axis(imagehandle,[0 max(xydata(:,1))*1.1 0 max(xydata(:,2))*1.1]);

%Run GUI
[xpos,ypos,tracks] = draw_maze_fcn(guihandle,fighandle,nexthandle,donehandle);