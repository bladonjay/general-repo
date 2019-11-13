function rotationGUI(tiffile)
    %# read image
    I = imread(tiffile);

    %# setup GUI
    hFig = figure('menu','none');
    hAx = axes('Parent',hFig);
    uicontrol('Parent',hFig, 'Style','slider', 'Value',0, 'Min',0,...
        'Max',360, 'SliderStep',[1 10]./360, ...
        'Position',[150 5 300 20], 'Callback',@slider_callback) 
    hTxt = uicontrol('Style','text', 'Position',[290 28 20 15], 'String','0');

    %# show image
    imshow(I, 'Parent',hAx)
    
    uicontrol('parent',hfig,'pushbutton','OK','position',[10 10 50 10],...
        'Callback',@getout)

    %# Callback function
    function slider_callback(hObj, eventdata)
        angle = round(get(hObj,'Value'));        %# get rotation angle in degrees
        imshow(imrotate(I,angle), 'Parent',hAx)  %# rotate image
        set(hTxt, 'String',num2str(angle))       %# update text
    end

    
end