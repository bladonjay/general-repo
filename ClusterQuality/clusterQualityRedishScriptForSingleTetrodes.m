% Script Implementing Redish's Cluster Quality Measures
% Implementing L_ratio and Isolation Distance from Schmitzer-Torbert et al.
% 2005 Nueroscience paper.

% Load data in from plexon, grab what you need.
% allChannelChannelIndex = NS15rec19091515001jakecutTetrode1ChannelFeatures(:,1)';
% unique(channelIndex);
% allChannelFeatureMatrix = NS15rec19091515001jakecutTetrode1ChannelFeatures(:,4:end)';
% allChannelUnitIndex = NS15rec19091515001jakecutTetrode1ChannelFeatures(:,2)';
% rawFeatureMatrix = allChannelFeatureMatrix(:,channelIndex==1);
% unitIndex = allChannelUnitIndex(:,channelIndex==1);

% Want all features, nothing that is labels.
% featureMatrix = channelFeaturesMatrix(:,4:end)';
% unitIndex = ChannelFeaturesNS15rec19091515001jakecut(:,2)';


%% Pattern for running things.
% rawFeatureMatrix = allChannelFeatureMatrix(:,allChannelChannelIndex==8);
% unitIndex = allChannelUnitIndex(:,allChannelChannelIndex==8);
% clusterQualityRedish
% isolationDistanceTetrode8 = isolationDistance;
% lRatioTetrode8 = lRatio;
% lTetrode8 = l;
% mahalanobisDistanceTetrode8 = mahalanobisDistance;
% featureMatrixTetrode8 = featureMatrix;
% unitIndexTetrode8 = unitIndex;

%% Calc mahalanobis distance
nRawFeatures = size(rawFeatureMatrix,1);
nClusters = length(unique(unitIndex))-1;
nSpikes = size(rawFeatureMatrix,2);
featureMatrix = zeros(0,nSpikes);
for i = 1:nRawFeatures
    if any(rawFeatureMatrix(i,:))
        featureMatrix(end+1,:) = rawFeatureMatrix(i,:);
    end
end
nFeatures = size(featureMatrix,1);

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

%% Calc measures of cluster quality
isolationDistance = nan(nClusters,1);
l = nan(nClusters,1);
lRatio = nan(nClusters,1);
for j = 1:nClusters
    notClusterJDistancesSorted = sort(mahalanobisDistance(unitIndex~=j,j),1,'ascend');
    l(j) = sum(chi2cdf(mahalanobisDistance(unitIndex~=j,j),nFeatures,'upper'));
    lRatio(j) = l(j)/sum(unitIndex==j);
    if sum(unitIndex==j) <= sum(unitIndex~=j)
        isolationDistance(j) = min(notClusterJDistancesSorted(sum(unitIndex==j)),...
            notClusterJDistancesSorted(round(0.1*sum(unitIndex~=j))));
    else
        isolationDistance(j) = notClusterJDistancesSorted(round(0.2*sum(unitIndex~=j)));
    end
end



%% Plotting
figure;
for j = 1:nClusters
    [ay,ax] = hist(mahalanobisDistance(unitIndex==j,j),10.^[0.1:0.1:5]);
    [ayNot,axNot] = hist(mahalanobisDistance(unitIndex~=j,j),10.^[0.1:0.1:5]);
    subplot(nClusters,1,j)
    plot(log10(ax),ay,'k')
    hold on;
    plot(log10(axNot),ayNot,'r')
end

%% Saving things
% isolationDistanceTetrode1 = isolationDistance;
% lRatioTetrode1 = lRatio;
% lTetrode1 = l;
% mahalanobisDistanceTetrode1 = mahalanobisDistance;
% featureMatrixTetrode1 = featureMatrix;
% unitIndexTetrode1 = unitIndex;
