function finalimage=rotationGUI2(picture,orig)
% picture is the image you want to rotate, orig is the image you want as
% your backghround reference, and finalimage is the rotated image.

% This function will take the canvas size of the original and plot 
% your image ontop of it- dont resize because pixels are your ground truth
% for image size

% updated from online rip to output a new image file at whatever image type
% the input is (rgb or one color)
%# read image
newimage=[];

% load your overlay, that will be the window size
[windowsize(1), windowsize(2)]=size(orig);

% make your image the same size
[overlay(1),overlay(2)]=size(picture);
%if windowsize(1)>overlay(1)
    

I = Grow2Square(picture);
% load your reference
I2 = Grow2Square(orig);
%# setup GUI
hFig = figure('menu','none');
hAx = axes('Parent',hFig);

% build the slider
uicontrol('Parent',hFig, 'Style','slider', 'Value',0, 'Min',-180,...
    'Max',180, 'SliderStep',[1 10]./360, ...
    'Position',[150 5 300 20], 'Callback',@slider_callback)


% put in the angle you have
hTxt = uicontrol('Style','text', 'Position',[290 28 20 15], 'String','0');

%# show initial image
%I2(:,:,1:2)=0; % set to red channel
imshow(I2, 'Parent',hAx)
hold on

% put an ok button so you can set angle, button calls uiresume
hWait=uicontrol('Parent',hFig,'Style','pushbutton','String','OK','Position',[10 10 50 20]...
    , 'callback','uiresume');
   

origsize=size(I);

% get the size
xhalf=origsize(1)/2; yhalf=origsize(2)/2;


% wait for the pushbutton to be pressed
uiwait

lastangle=str2double(get(hTxt, 'String'));
finalimage=imrotate(I,lastangle);
% resize

xhalf=fix(origsize(1)/2); yhalf=fix(origsize(2)/2);
lastsize=size(finalimage);
% crop at center
xindices1=[fix(lastsize(1)/2)-xhalf fix(lastsize(1)/2)+xhalf];
yindices1=[fix(lastsize(2)/2)-yhalf fix(lastsize(2)/2)+yhalf];
finalimage=finalimage(xindices1(1):xindices1(2),yindices1(1):yindices1(2),:);
close;



%# Callback function
    function slider_callback(hObj, eventdata)
        angle = round(get(hObj,'Value'));        %# get rotation angle in degrees
        
        % rotate image
        
        tempimage=imrotate(I,angle);
        %tempimage(:,:,2:3)=0;
        %white=max(max(max(tempimage)));
        %tempimage(:,:,2:3)=white;  % set to different color;
        newsize=size(tempimage);
        % crop at center, usually the new size isnt smaller than the old
        % but you never know...
        xindices=[ceil(newsize(1)/2)-xhalf floor(newsize(1)/2)+xhalf];
        yindices=[ceil(newsize(2)/2)-yhalf floor(newsize(2)/2)+yhalf];
        tempimage=tempimage(xindices(1):xindices(2),yindices(1):yindices(2),:);
        %I2(:,:,1)=white;
        imshow(I2,'Parent', hAx);
        J=imshow(tempimage, 'Parent',hAx);  %# rotate image
        % set alphadata to let you see underlying image
        set(J,'AlphaData',tempimage(:,:,3));
        set(hTxt, 'String',num2str(angle))       %# update text
    end


end