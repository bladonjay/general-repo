%% Matching up clusters from cluster quality work with the corresponding 
%  neurons in our dataset.

%    To create a file list:
%       DirData = dir('ClusterQualityData\');
%       fileList = {DirData(3:end).name}';

% ClusterQualitySubPaperMainDataset = clusterQualityMetrics(fileList);
% load('D:\AxisCellPaper\Analyses\WaveformClusteringAnalysis\spikeFormData_complete.mat');

ABC = spikeFormData;
XYZ = ClusterQualitySubPaperMainDataset;

for i = 1:length(ABC.channel)
    ABC.intChannel(i) = str2double(ABC.channel{i}(5:6));
    ABC.intNeuron(i) = double(ABC.channel{i}(7)-96);
end
clusterLabels = [XYZ.rat;XYZ.rec;XYZ.channel;XYZ.neuronNumber]';
neuronsAnalyzedLabels = [ABC.rat;ABC.rec;ABC.intChannel;ABC.intNeuron]';

indexInClusterStruct = zeros(542,1);
for i = 1:542;
    if sum(ismember(clusterLabels,neuronsAnalyzedLabels(i,:),'rows')) > 1
        break
    elseif any(ismember(clusterLabels,neuronsAnalyzedLabels(i,:),'rows'))
        indexInClusterStruct(i) = find(ismember(clusterLabels,neuronsAnalyzedLabels(i,:),'rows'));
    end
end

find(indexInClusterStruct == 0)

% Then you have to go in by hand and fill in any unlabeled due to
% discrepancies in the wire of the tetrode used to label the tetrode (file 
% naming errors cause this).


