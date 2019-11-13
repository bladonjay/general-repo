function CheckPos(X,Y,ts)
% I think this just plots a figure in which you can grab the ts for the
% datapoint.

plot( X,Y,'.','color','white','markersize',.1)


dcm = datacursormode(gcf);
datacursormode on;
set(dcm, 'updatefcn', @myfunction)




function output_txt = myfunction( obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
% event_obj

dataIndex = get(event_obj,'DataIndex');
pos = get(event_obj,'Position');

output_txt = {[ 'X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};

    output_txt{end+1} = ['ts =' num2str(ts(dataIndex))];



% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end
end
end