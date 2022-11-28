function [ ClassicSpatialMeasures ] = classicSpatialMeasures(RecStruct, NeuronStruct)
%   Calcs spatial info, selectivity, sparsity, coherence
%
%
%   Written by Jake Olson, October 2015
%   Last updated by Jake Olson, Oct 2020.

% Coherence, Selectivity, Info per Spike
MIN_OCC_SPATIAL_STATS = 1/10;
SAMPLE_RATE = 60; % same for all our recs
ClassicSpatialMeasures.minOccForStatInclusionSecs = MIN_OCC_SPATIAL_STATS;

PIXELSPERCM = 3.5;
BINSIZECM = 2;
minOccInSamples = SAMPLE_RATE * MIN_OCC_SPATIAL_STATS;  % 1/10 of a second.

% max bin values can change w/out being considered a jump and removed
JUMP_MAX = 35;

for iNeu = 1:length(NeuronStruct.neuronName)
    recIndex = NeuronStruct.recStructIndex(iNeu);
    
    % get spike and tracking data, remove bad tracking
    if isnan(RecStruct.sessionTimeStamps{recIndex})
        startEpochSample = 1;
        endEpochSample = length(RecStruct.Behavior.posXY{recIndex});
    else
        startEpochSample = RecStruct.sessionTimeStamps{recIndex}(3);
        endEpochSample = RecStruct.sessionTimeStamps{recIndex}(4);
    end
    tracking = RecStruct.Behavior.posXY{recIndex}(startEpochSample:endEpochSample,:);
    spikes = NeuronStruct.sampleNSpikes{iNeu}(startEpochSample:endEpochSample);
    
    xs = round(tracking(:,1) / (BINSIZECM * PIXELSPERCM))+1; %adding 1 avoids 0s, which can cause problems.
    ys = round(tracking(:,2) / (BINSIZECM * PIXELSPERCM))+1;
    %take out rows with lost tracking
    badSpotsLost = isnan(xs) | isnan(ys);
    % remove jumps from glare from plotting.
    badSpotsJumps = [0; abs(diff(xs)) > JUMP_MAX | abs(diff(ys)) > JUMP_MAX ];
    anyBadSpots = badSpotsLost | badSpotsJumps;
    xsClean = xs(~anyBadSpots);
    ysClean = ys(~anyBadSpots);
    xMax = max(xsClean);
    yMax = max(ysClean);
    
    % 2D rate map
    [frMap, occMap] = mapVar2D(spikes(~anyBadSpots), xsClean, ysClean, [xMax,yMax]);
    frMap = frMap * SAMPLE_RATE;
    
    % Only include places that were occupied at least x times.
    nOccs = sum(sum(occMap(occMap >= minOccInSamples)));
    pBin = occMap./nOccs;
    pBin(occMap < minOccInSamples) = 0;
    goodBins = (pBin ~= 0);
    
    meanFR = sum(sum(frMap(goodBins).*pBin(goodBins)));
    maxFR =  max(max(frMap(goodBins)));
    
    ClassicSpatialMeasures.meanFR(iNeu) = meanFR;
    ClassicSpatialMeasures.maxFR(iNeu) = maxFR;
            
    % if no firing on track:
    if meanFR == 0
        ClassicSpatialMeasures.maxFR(iNeu) = 0;
        ClassicSpatialMeasures.spatialInfoPerSpikeInBits(iNeu) = NaN;
        ClassicSpatialMeasures.spatialInfoInBits(iNeu) = NaN;
        ClassicSpatialMeasures.sparsity(iNeu) = NaN;
        ClassicSpatialMeasures.selectivity(iNeu) = NaN;
        ClassicSpatialMeasures.spatialCoherence(iNeu) = NaN;
        ClassicSpatialMeasures.spatialCoherencePVal(iNeu) = NaN;
    else

        %spatial information per spike - 1993 Skaggs et al. NIPS
        tmpSpInfo = (pBin(goodBins).*...
            frMap(goodBins)/meanFR).*...
            log2(frMap(goodBins)/meanFR);
        ClassicSpatialMeasures.spatialInfoPerSpikeInBits(iNeu) = ...
            sum(sum(tmpSpInfo(~isnan(tmpSpInfo))));
        
        %spatial information per second- 1993 Skaggs et al. NIPS
        tmpSpInfo = (pBin(goodBins).*...
            frMap(goodBins)).*...
            log2(frMap(goodBins)/meanFR);
        ClassicSpatialMeasures.spatialInfoInBits(iNeu) = ...
            sum(sum(tmpSpInfo(~isnan(tmpSpInfo))));
        
        % spatial sparsity - 1996 Skaggs  et al. hippocampus
        ClassicSpatialMeasures.sparsity(iNeu) = ...
            nansum(nansum(pBin.*frMap))^2/...
            nansum(nansum(pBin.*frMap.^2));
        
        % spatial selectivity - 1996 Skaggs  et al. hippocampus
        ClassicSpatialMeasures.selectivity(iNeu) = maxFR/meanFR;
        
        % Spatial Coherence - from Andy - a la Kubie Muller Bostock 1990
        countBin = 0;
        neighborFR = [];
        actualFR = [];
        for x = 2:size(frMap,1)-1
            for y = 2:size(frMap,2)-1
                if occMap(x,y) ~= 0
                    countBin = countBin + 1;
                    neighborhood = frMap(x-1:x+1,y-1:y+1);
                    neighborFR(countBin) = nanmean(neighborhood([1:4,6:9]));
                    actualFR(countBin) = frMap(x,y);
                end
            end
        end
        [ClassicSpatialMeasures.spatialCoherence(iNeu)] = ...
            corr(actualFR',neighborFR','rows','complete','type','Spearman','tail','right');
    end
    
end





