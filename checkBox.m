function checked = checkBox(markers,name)
% function [checked]=checkBox(markers,name);
% markers is cell array of words to go in each box
% name is a string
% checked will be a set of indices corresponding to those that were checked

if ~exist('name','var')
    name='choose one or some';
end
checked=[];
screen = get(0,'ScreenSize');
% Create figure
boxsizey=max([30*length(markers)+30, 400]);
boxsizex=max([150 max(cellfun(@length, markers))*15]);

h.f = figure('units','pixels','position',[200,(screen(4)-boxsizey-400),boxsizex,boxsizey],...
    'toolbar','none','menu','none','Name',name);


%Create yes/no checkboxes
for i = 1:length(markers)
    h.c(i) = uicontrol('style','checkbox','units','pixels',...
        'position',[10,(i-1)*30 + 30,boxsizex,15],'string',markers{i});
end

% vals = get(h.c(:),'Value');
% checked = find([vals{:}]);
%
% if isempty(checked)
%    checked = [];
% end



% Create OK pushbutton
h.p = uicontrol('style','pushbutton','units','pixels',...
    'position',[40,5,70,20],'string','OK', ...
    'Callback',@p_call);


%fprintf('asdf');
uiwait(h.f);

    function p_call(c,varargin)
        vals = get(h.c(:),'Value');
        
        if iscell(vals)
            checked = find([vals{:}]);
        else
            checked=vals;
        end
        
        if isempty(checked)
            checked = [];
        end
        close
        
    end


end

%
% function p_call(c,varargin)
% vals = get(h.c(:),'Value');
%
% checked = find([vals{:}]);
%
% if isempty(checked)
%    checked = [];
% end
% close
%
% end
