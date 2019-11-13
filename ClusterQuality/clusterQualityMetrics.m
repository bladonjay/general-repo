function [ ClusterQualityStruct ] = clusterQualityMetrics( fileList )
%CLUSTERQUALITYMETRICS Implementing Redish's Cluster Quality Measures
%   Implementing L_ratio and Isolation Distance from Schmitzer-Torbert et al.
%   2005 Nueroscience paper.
%
%    [ClusterQualityStruct] =  CLUSTERQUALITYMETRICS( fileList ) takes one
%       argument, a cell array with a list of file names to run the cluster
%       quality analyses on. If a simpler one tetrode analysis is wanted,
%       use clusterQualityRedishScriptForSingleTetrodes. It is the same
%       code without the loop and function wrapper.
%
%    OUTPUT: 1 Struct with 11 fields.
%    ClusterQualityStruct is a struct with 11 fields. Each field is a vector
%       of size N where N is the total number of neurons processed. The
%       fields are:
%           rat - rat number - all NS
%           rec - recording number
%           channel - channel number of first open wire of the tetrode
%           neuronNumber - index of the neuron in the tetrode. Numbers
%               correspond to the matching letter in offlineSorter. with
%               the first 4 fields, the exact neuron can be identified and
%               matched to other analyses.
%           l - Redish cluster quality measure. see 2005 paper. Correlates
%               well with type 2 errors (improperly excluded spikes).
%           lRatio - l value normalized by number of spikes in the cluster.
%               This is a better measure than L since it accounts for
%               expecting more noise/misses with larger clusters.
%           isolationDistance - A Harris and Buzsaki (2001) measure, it
%               finds the Mahalanobis distance of the Nth noise spike from
%               the cluster center as evaluated by Mahalanobis distance,
%               where N is the number of spikes in the cluster being
%               evaluated. Since this is undefined when the cluster has
%               more spikes than the noise population, and starts to lose
%               intuitive meaning approaching that state, I have adapted
%               isolationDistance to use either the distance of Nth noise
%               spike or the distance of the noise spike 20% into the noise
%               spike distance distribution (the 200th closest if 1k noise 
%               spikes), whichever is closer. This reduces the score for
%               large clusters and I believe it is more representative of 
%               the cluster quality.
%           nFeatures - Mahalanobis distance scales with dimensionality,
%               and our feature number can change, especially if a channel
%               on the tetrode was closed due to noise. This is the number
%               of features used in the calculation.
%           nNeurons - The number of neurons on this tetrode. Mostly useful
%               only as a check of veracity of the tetrode files to the
%               actual data.
%           nSpikes - The number of spikes for the cluster. Same purpose
%               as nNeurons.
%           adaptedIsoDistFlag - 1 if the neuron used the adapted
%               IsolationDistance for many spike neurons detailed above, 0 
%               if not.
%
%    To create a file list:
%       DirData = dir('ClusterQualityData\');
%       fileList = {DirData(3:end).name}';
%
%    Note: What if there are no noise spikes? then put in NaNs.
%
%    Non-built-in functions called:
%
%   Written by Jake Olson, July 2016
%   Last updated by Jake Olson, July 2016.




% Get number of files to process. All should be tab delimited text files
% with a header row. Unit index is the first nontext feature column.
nTetrodes = length(fileList);
iNeuron = 0;
for iTetrode = 1:nTetrodes
    % for each tetrode, import waveform feature data and cluster label data
    waveformFeatureMat = importdata(fileList{iTetrode},'\t',1);
    rawFeatureMatrix = waveformFeatureMat.data(:,2:end)';
    unitIndex = waveformFeatureMat.data(:,1)';
    
    % Define variables based on matrix sizes.
    nRawFeatures = size(rawFeatureMatrix,1);
    nClusters = sum(unique(unitIndex)>0);
    nSpikes = size(rawFeatureMatrix,2);
    
    % Remove features that are all 0 (no data, due to closed recording
    % channel.
    featureMatrix = zeros(0,nSpikes);
    for i = 1:nRawFeatures
        if any(rawFeatureMatrix(i,:))
            featureMatrix(end+1,:) = rawFeatureMatrix(i,:);
        end
    end
    nFeatures = size(featureMatrix,1);
    
    % Calc mahalanobis distance
    mahalanobisDistance = nan(nSpikes,nClusters);
    for j = 1:nClusters
        mahalanobisDistance(:,j) = mahal(featureMatrix',featureMatrix(:,unitIndex==j)');
        % same as built in matlab line above - but slower.
        % meanFeaturesClusterC = mean(featureMatrix(:,unitIndex==j),2);
        % covMatClusterC = cov(featureMatrix(:,unitIndex==j)');
        % for i = 1:nSpikes
        %     mahalanobisDistance(i,j) = (featureMatrix(:,i)-meanFeaturesClusterC)'/covMatClusterC*(featureMatrix(:,i)-meanFeaturesClusterC);
        % end
    end
    
    % Calc measures of cluster quality
    isolationDistance = nan(nClusters,1);
    l = nan(nClusters,1);
    lRatio = nan(nClusters,1);
    for j = 1:nClusters
        if sum(unitIndex~=j) == 0 % one cluster, no noise spikes
            iNeuron = iNeuron+1;
            ClusterQualityStruct.l(iNeuron) = NaN;
            ClusterQualityStruct.lRatio(iNeuron) = NaN;
            ClusterQualityStruct.isolationDistance(iNeuron) = NaN;
            ClusterQualityStruct.rat(iNeuron) = str2double(fileList{iTetrode}(3:4));
            ClusterQualityStruct.rec(iNeuron) = str2double(fileList{iTetrode}(9:10));
            ClusterQualityStruct.channel(iNeuron) = str2double(fileList{iTetrode}(12:13));
            ClusterQualityStruct.neuronNumber(iNeuron) = j;
            ClusterQualityStruct.nFeatures = nFeatures;
            ClusterQualityStruct.nNeurons = nClusters;
            ClusterQualityStruct.nSpikes = sum(unitIndex==j);
            ClusterQualityStruct.adaptedIsoDistFlag = 0;
            break
        end
        notClusterJDistancesSorted = sort(mahalanobisDistance(unitIndex~=j,j),1,'ascend');
        l(j) = sum(chi2cdf(mahalanobisDistance(unitIndex~=j,j),nFeatures,'upper'));
        lRatio(j) = l(j)/sum(unitIndex==j);
        if 5*sum(unitIndex==j) <= sum(unitIndex~=j)
            isolationDistance(j) = notClusterJDistancesSorted(sum(unitIndex==j));
            adaptedIsoDistFlag = 0;
        else
            isolationDistance(j) = notClusterJDistancesSorted(round(0.2*sum(unitIndex~=j)));
            adaptedIsoDistFlag = 1;
        end
        
        % Save results
        iNeuron = iNeuron+1;
        ClusterQualityStruct.l(iNeuron) = l(j);
        ClusterQualityStruct.lRatio(iNeuron) = lRatio(j);
        ClusterQualityStruct.isolationDistance(iNeuron) = isolationDistance(j);
        ClusterQualityStruct.rat(iNeuron) = str2double(fileList{iTetrode}(3:4));
        ClusterQualityStruct.rec(iNeuron) = str2double(fileList{iTetrode}(9:10));
        ClusterQualityStruct.channel(iNeuron) = str2double(fileList{iTetrode}(12:13));
        ClusterQualityStruct.neuronNumber(iNeuron) = j;
        ClusterQualityStruct.nFeatures = nFeatures;
        ClusterQualityStruct.nNeurons = nClusters;
        ClusterQualityStruct.nSpikes = sum(unitIndex==j);
        ClusterQualityStruct.adaptedIsoDistFlag = adaptedIsoDistFlag;
    end
end

%% Plotting
% figure;
% for j = 1:nClusters
%     [ay,ax] = hist(mahalanobisDistance(unitIndex==j,j),10.^[0.1:0.1:5]);
%     [ayNot,axNot] = hist(mahalanobisDistance(unitIndex~=j,j),10.^[0.1:0.1:5]);
%     subplot(nClusters,1,j)
%     plot(log10(ax),ay,'k')
%     hold on;
%     plot(log10(axNot),ayNot,'r')
% end
%

end

