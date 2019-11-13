% ginputSmooth: continuously store mouse position while left mouse button is held down
%
% Returns list of x and y coordinates in axes units
% Ideal for use in MClust cutter utility (e.g. DrawConvexHull) instead of ginput
%
% Written by BRK 2015. Please distribute and modify freely.

function [X,Y] = ginputSmooth()
% need zoom off to use windowbuttonfcns
zoom off
% initialize coordinate vectors and handles structure
X = [];
Y = [];
finalData = [];
dynamicData = guihandles(gcf);
dynamicData.x = [];
dynamicData.y = [];
guidata(gcf,dynamicData);
figName = get(gcf,'name');
% set button functions
set(gcf,'windowbuttondownfcn',@click);
    function click(~,~)
        % clicking activates motion and release functions
        set(gcf,'windowbuttonmotionfcn',@move);
        set(gcf,'windowbuttonupfcn',@release);
        function move(~,~)
            % store mouse coordinates as it moves
            currentPos = get(gca,'CurrentPoint');
            currentX = currentPos(1,1);
            currentY = currentPos(1,2);
            dynamicData = guidata(gcbo);
            dynamicData.x = [dynamicData.x; currentX];
            dynamicData.y = [dynamicData.y; currentY];
            guidata(gcbo,dynamicData);
        end
        function release(~,~)
            % reset button fcns and get output
            set(gcf,'windowbuttonmotionfcn','');
            set(gcf,'windowbuttondownfcn','');
            finalData = guidata(gcbo);
            X = finalData.x;
            Y = finalData.y;
            % need to alter fig property to trigger resume of waitfor
            set(gcf,'name','cut');
        end
    end
% wait for fig property change triggered by button release
waitfor(gcf,'name','cut')
% reset fig name
set(gcf,'name',figName)
end

