
function [fixed,brightness]=brightengui(orig,method)


if (~exist('method','var') || ~ischar(method))
    method = nan;
end

if strcmpi(method,'multiply')
    scaler=1;
elseif strcmpi(method,'bottom')
    scaler=2;
elseif strcmpi(method,'top')
    scaler=3;
elseif strcmpi(method,'add')
    scaler=4;
else
    scaler=1;
end


tempimage=[];
%if method=='multiply'

% build figure with original
hFig = figure('menu','none');
hAx = axes('Parent',hFig);

%   

% build the slider
uicontrol('Parent',hFig, 'Style','slider', 'Value', 20, 'Min',20,...
    'Max',100, 'SliderStep',[.01, .90], ...
    'Position',[150 5 300 20], 'Callback',@slider_callback)

% put in the brightness you have
hTxt = uicontrol('Style','text', 'Position',[290 28 20 15], 'String','0');


image(orig, 'Parent',hAx)
axis off;

hWait=uicontrol('Parent',hFig,'Style','pushbutton','String','OK','Position',[10 10 50 20]...
    , 'callback','uiresume');

uiwait
fixed=tempimage;
% dont really need to plot it out
%image(fixed, 'Parent',hAx)
close;

%# Callback function
    function slider_callback(hObj, eventdata)
        brightness = get(hObj,'Value');        %# get rotation angle in degrees
        
        % fix image using method

        switch scaler
            case 1 % multiply
                tempimage=orig*(brightness/20);
            case 2 % squeeze bottom
                oldx=min(linearize(orig)) : max(linearize(orig));
                newx=min(linearize(orig))+brightness : max(linearize(orig));
                [~,~,c]=histcounts(linearize(orig),length(newx));
                newpx=newx(c);
                 tempimage=reshape(newpx,size(orig,1),size(orig,2),size(orig,3));
            case 3 % squeeze top
                oldx=min(linearize(orig)) : max(linearize(orig));
                newx=min(linearize(orig)) : max(linearize(orig)-brightness);
                [~,~,c]=histcounts(linearize(orig),length(newx));
                newpx=newx(c);
                tempimage=reshape(newpx,size(orig,1),size(orig,2),size(orig,3));
            case 4 % add to everyone
                tempimage=orig+brightness;
        end

        %I2(:,:,1)=white;
        image(tempimage,'Parent', hAx);
        axis off;
        % set alphadata to let you see underlying image

        set(hTxt, 'String',num2str(brightness))       %# update text
    end

end


