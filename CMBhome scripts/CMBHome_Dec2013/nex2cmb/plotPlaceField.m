
%get epochs
%%
doorIDX = cellfun(@any,strfind(root.event(:,1),'Door Open'));
rewIDX = cellfun(@any,strfind(root.event(:,1),'Reward'));
doorCloseIDX = cellfun(@any,strfind(root.event(:,1),'Door Closed'));

doorTS = cell2mat(root.event(doorIDX,2));
reward = cell2mat(root.event(rewIDX,2));
doorclosed = cell2mat(root.event(doorCloseIDX,2));

%%

if length(doorTS) == length(doorclosed)
   epochs = [ doorTS doorclosed];
   cor = cellfun(@(a) any(reward >=a(1) & reward <=a(2)),num2cell(epochs,2));
else
    keyboard
end

a = ones(length(doorTS),1);
b = 2*ones(length(reward),1);
c = 3*ones(length(doorclosed),1);
match = [a; b;c];
[ts,idx] = sort([doorTS;reward;doorclosed]);
match = match(idx);


i = find(diff([0;match]) == 0);

for j = 1:length(i)
   ev = match(i(j));
    switch ev
        case 1
            doorTS(ismember(doorTS,ts(i(j)))) = [];
        case 2
             reward(ismember(reward,ts(i(j)))) = [];
        case 3
             doorclosed(ismember(doorclosed,ts(i(j)))) = [];
            
    end
end

%%
root = root.FixPos;

rang_x = [min(root.b_x) max(root.b_x)];
rang_y = [min(root.b_y) max(root.b_y)];

binX = rang_x(1):12:rang_x(2);
binY = rang_y(1):12:rang_y(2);

for i = 1:size(root.cells,1)
    figure
   spk_x = cell2mat(root.spk_x(root.cells(i,:)));
   spk_y = cell2mat(root.spk_y(root.cells(i,:)));
    
   
   heatmap = histcn([spk_x spk_y],binX,binY);
    imagesc(heatmap)
    colorbar
end


%%


binnedPop=populationMatrix(root,3,3,60,reward,0);