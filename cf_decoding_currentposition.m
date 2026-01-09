function cf_decoding_currentposition(dirs, params, allindex, metadata, posType, cellTypes)
%cf_decoding_currentposition
%
%ALP 12/28/2022

%%% group vectors
isGamma = strcmp(metadata.Groups, 'gamma'); 
isRandom = strcmp(metadata.Groups, 'random');
isPre = metadata.FlickerDay < 5;
isPost = metadata.FlickerDay > 5; 

%%% params
params.train_perc = 80; 
params.test_perc = 100-params.train_perc; 
params.k_subsets = 100/params.test_perc; 
params.minTrials = params.k_subsets;
if strcmp(posType, 'full')
    params.posEdges = 0:6:360;
    params.posBins = 6;
    dType = '360';
    nDeg = 360;
    posName = 'theta';
    params.RZ = [54 72; 234 252];
else
    params.posEdges = -81:3:99;
    params.posBins = 3;
    dType = '180';
    nDeg = 180; 
    posName = 'theta_d2r'; 
    params.RZ = [0 18];
end
params.cellTypes = cellTypes; 
params.decodingBins = 0.2; %in s
params.speedThreshold = 1; 

%%% directories 
savedatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\decoding_currentposition\';
trialdatadir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\positionInfo\'; 
ratemapdir = '\\ad.gatech.edu\bme\labs\singer\Abby\code\chronicflicker-ephys-prospectivecoding\results\placecells\';

%%% filenames
trialdatafilename = ['trialInfo_', dType, '_A'];
trialspikesfilename = ['trialSpikes_', dType, '_A'];
savefilename = ['decoding_currpos_', dType, '_', params.cellTypes, '_test', num2str(params.test_perc)]; 

%%% set things as needed
dayindex = unique(allindex(:,1:2), 'rows'); 

%%% loop over all days
for d = 1:size(dayindex,1)
    %%% load trial data
    load([trialdatadir, trialdatafilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
    load([trialdatadir, trialspikesfilename, num2str(dayindex(d,1)), '_', num2str(dayindex(d,2)), '.mat'])
    
    %%% get trials to include
    inclTrials = [trialData.engaged] & [trialData.fullTrial] & [trialData.rewarded]; 
    goodTrials = trialData(inclTrials);
    [goodTrials.spikeTimes] = deal(trialSpikes.data(inclTrials).spikeTimes);
    [goodTrials.spikePosInds] = deal(trialSpikes.data(inclTrials).spikePosInds);
    [goodTrials.spikeIDs] = deal(trialSpikes.data(inclTrials).IDs); 
    
    if strcmp(params.cellTypes, 'all')
        inclCells = trialSpikes.info.clusterID;
    elseif strcmp(params.cellTypes, 'PC')
        load([ratemapdir, 'spatialmaps_all_500_full_', num2str(dayindex(d,2)), '.mat'])
        inclPC = ratemaps.includeCell;
        inclCells = trialSpikes.info.clusterID(inclPC);
    end
    
    %%% check if enough trials to include this day
    if length(goodTrials) < params.minTrials
        disp(['skipping day ', num2str(dayindex(d,2)), ' bc not enough trials'])
        groupData(d).nTrials = length(goodTrials);
        continue
    end
    
    %%% get testing and training indices
    kInds = crossvalind('kfold', length(goodTrials), params.k_subsets);
    
    for iK = 1:params.k_subsets
        %trainingData
        allIDs = inclCells;
        nCells = length(allIDs);
        isTrainTrial = kInds ~= iK;
        trainingData = cf_gettrainingdata_currentposition(goodTrials(isTrainTrial), allIDs, nCells, params.speedThreshold, params.posEdges, posName);
        
        %testingData
        isTestTrial = kInds == iK;
        testingData = cf_gettestingdata_currentposition(goodTrials(isTestTrial), allIDs, nCells, params.speedThreshold, params.decodingBins, posName);
        
        %decode each testing trial
        decodedData(iK) = cf_getdecodedposition_currentposition(trainingData, testingData, nCells, params.decodingBins, params.posEdges, posName);
    
        %append across kFolds
        mnEstPos(:,:,iK) = decodedData(iK).mnEstPos;
        binError(iK,:) = decodedData(iK).binError; 
    end
    
    %%% get mean and diagonals
    dayMean = mean(mnEstPos,3,'omitnan');
    diagVals = NaN(3, length(params.posEdges(1:end-1))); 
    diagVals(1,:) = [NaN diag(dayMean,1)'];
    diagVals(2,:) = diag(dayMean);
    diagVals(3,:) = [diag(dayMean, -1)', NaN];
    diagErr = mean(mean(diagVals,1, 'omitnan'), 'omitnan'); 
    
    figure
    hold on
    sgtitle({['curr position decoding ', num2str(dayindex(d,2)), ' ', dType, ' celltypes ', params.cellTypes, ' diagErr ', num2str(diagErr)], [' nTrials ', num2str(length(goodTrials)), ' nCells ', num2str(nCells)]})
    subplot(2,1,1)
    hold on
    if ~isnan(diagErr)
        clims = [0 0.75*max(max(dayMean))];
        imagesc(params.posEdges, params.posEdges, dayMean, clims)
    else
        imagesc(params.posEdges, params.posEdges, dayMean)
    end
    RZxplot = sort(reshape(params.RZ, [numel(params.RZ),1]));
    RZxplot = repmat(RZxplot, [1, length(params.posEdges(1:end-1))]);
    RZyplot = repmat(params.posEdges(1:end-1), [numel(params.RZ),1]); 
    plot(RZxplot', RZyplot', 'r--')
    xlim([min(params.posEdges) max(params.posEdges)])
    ylim([min(params.posEdges) max(params.posEdges)])
    xlabel('actual position (deg)')
    ylabel('estimated position (deg)')
    colorbar
    subplot(2,1,2)
    hold on
    plot(params.posEdges(1:end-1), mean(binError,1, 'omitnan'))
    plot(RZxplot', RZyplot', 'r--')
    xlim([min(params.posEdges) max(params.posEdges)])
    ylim([-nDeg/4 nDeg/4])
    xlabel('actual position (deg)')
    ylabel('decoding error, + in front')
    savefigfilename = ['fig_decoding_currpos_', dType, '_', params.cellTypes, '_test', num2str(params.test_perc) '_', num2str(dayindex(d,2))];
    savefigALP(savedatadir, savefigfilename, 'filetype', 'pdf')

    decodingResult.mnEstPos = dayMean;
    decodingResult.decodedData = decodedData; 
    decodingResult.binError = mean(binError,1, 'omitnan'); 
    decodingResult.diagonalError = diagErr; 
    
    info = []; 
    info = addhelpfulinfotostruct(info); 
    save([savedatadir, savefilename, '_', num2str(dayindex(d,2)), '.mat'], 'decodingResult', 'info')
    
    groupData(d).mnEstPos = dayMean;
    groupData(d).binError = mean(binError,1,'omitnan'); 
    groupData(d).avgBinError = mean(groupData(d).binError, 'omitnan');
    groupData(d).diagonalError = diagErr; 
    groupData(d).nTrials = length(goodTrials);
    groupData(d).nCells = nCells;
    
    clear decodingResult goodTrials decodedData trainingData testingData mnEstPos binError diagErr
end

%%% get days to include
inclDay = zeros(size(dayindex,1),1);
for d = 1:size(dayindex,1)
    if ~isempty(groupData(d).diagonalError)
        if groupData(d).diagonalError > 0.02
            inclDay(d) = 1;
        end
    end
    groupData(d).include = inclDay(d);
    
end

%%% get distribution of digonal values for inclusion
figure
hold on
allDiagErr = [groupData.diagonalError]; 
allDiagErr = allDiagErr.*(180/3);
histogram(allDiagErr, 0:0.1:10)

cellPosEst = {groupData.mnEstPos};
gammaPosEst = cat(3,cellPosEst{isGamma&isPost&inclDay}); 
randomPosEst = cat(3,cellPosEst{isRandom&isPost&inclDay}); 

figure('Position', [440 613 644 185])
hold on
title('probability/chance - post flicker sess - all sess')
subplot(1,2,1)
hold on
imagesc(params.posEdges, params.posEdges, mean(randomPosEst,3,'omitnan').*(nDeg/params.posBins), [0.25 3])
plot(RZxplot', RZyplot', 'w--')
xlim([min(params.posEdges) max(params.posEdges)])
ylim([min(params.posEdges) max(params.posEdges)])
xlabel('actual position (deg)')
ylabel('estimated position (deg)')
colorbar
title('random')
subplot(1,2,2)
hold on
imagesc(params.posEdges, params.posEdges, mean(gammaPosEst,3,'omitnan').*(nDeg/params.posBins), [0.25 3])
plot(RZxplot', RZyplot', 'w--')
xlim([min(params.posEdges) max(params.posEdges)])
ylim([min(params.posEdges) max(params.posEdges)])
xlabel('actual position (deg)')
ylabel('estimated position (deg)')
title('gamma')
colorbar
makefigurepretty(gcf)
savefigfilename = ['GROUP_decoding_currpos_', dType, '_', params.cellTypes, '_test', num2str(params.test_perc)];
savefigALP(savedatadir, savefigfilename)


gammaBinError = {groupData(isGamma&isPost&inclDay).binError};
gammaBinError = cell2mat(gammaBinError');
randomBinError = {groupData(isRandom&isPost&inclDay).binError};
randomBinError = cell2mat(randomBinError');

gammaMnError = mean(gammaBinError,1,'omitnan');
randomMnError = mean(randomBinError,1,'omitnan');
stdeGammaError = std(gammaBinError, [], 1, 'omitnan')./sqrt(sum(~isnan(gammaBinError(:,1))));
stdeRandomError = std(randomBinError, [], 1, 'omitnan')./sqrt(sum(~isnan(randomBinError(:,1))));

figure
hold on
shadedErrorBar(params.posEdges(1:end-1), gammaMnError, stdeGammaError, {'Color', params.colors.gamma.post},1)
shadedErrorBar(params.posEdges(1:end-1), randomMnError, stdeRandomError, {'Color', params.colors.random.post},1)
plot(RZxplot', RZyplot', 'r--')
xlim([min(params.posEdges) max(params.posEdges)])
ylim([-nDeg/8 nDeg/8])
xlabel('position (deg)')
ylabel('decoding error (deg)')
title('decoding error - all sess - post flicker')
makefigurepretty(gcf)
savefigfilename = ['GROUP_decodingError_currpos_', dType, '_', params.cellTypes, '_test', num2str(params.test_perc)];
savefigALP(savedatadir, savefigfilename)

savefilename = ['Group_decoding_currpos_', dType, '_', params.cellTypes, '_test', num2str(params.test_perc)]; 
save([savedatadir, savefilename, '.mat'], 'groupData', 'info')




end

