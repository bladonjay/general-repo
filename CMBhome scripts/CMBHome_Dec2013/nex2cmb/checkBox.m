function checked = checkBox(markers)
% Create figure


h.f = figure('units','pixels','position',[200,200,100,30*length(markers)+30],...
    'toolbar','none','menu','none');


for i = 1:length(markers)
    % Create yes/no checkboxes
    h.c(i) = uicontrol('style','checkbox','units','pixels',...
        'position',[10,(i-1)*30 + 30,100,15],'string',markers{i});
 
    
    
end
% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[40,5,70,20],'string','OK',...
    'callback',@p_call);
% Pushbutton callback
uiwait(h.f)
    function p_call(varargin)
        vals = get(h.c,'Value');
       
        checked = find([vals{:}]);
        
        if isempty(checked)
            checked = [];
        end
        close(h.f)
    end
end