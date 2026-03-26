function [decodedEstPosAll, decodedErrorAll] = getDecodingResultsByDay_linearJLK(dirs,uniqSess,params)
%adapted and integrated cf_RunAnalyses-->cf_decoding_currentposition-->...
% cf_gettrainingdata_currentposition.m/cf_gettestingdata_currentposition.m/...
% cf_getdecodedposition_currentposition.m/cf_getspatialprob
%code follows similar structure as getBehaviorROC_JLK.m,
%getFiringRateResultsByDay_linearJLK, and
%getSpeedandLickRateResultsByDay_linearJLK

%% perform decoding with day data %%

for id = 2%:length(params.decoding.decID)

    % % % %output info
    % % % decResultsfname = fullfile(dirs.saveoutputstructs, 'Data\Neural\Decoding', [params.decoding.decID{id} 'Results.mat']);
    % % 
    % % %initialize output struct
    % % decodedEstPosAll = [];
    % % decodedActualPosAll = [];
    % % decodedErrorAll = [];

    %loop through unique sessions
    for ss = 1%:51%size(uniqSess,1)

        %load day data
        dayDatafname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\JK' num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' params.decoding.decID{id} 'Data.mat']);
        decResultsDir = fullfile([dirs.saveoutputstructs, 'Data\Neural\Decoding\JK' num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2))]);
        decResultsfname = fullfile([dirs.saveoutputstructs, 'Data\Neural\Decoding\JK' num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' params.decoding.decID{id} 'DecData.mat']);

        if isfile(dayDatafname) && ~isfile(decResultsfname)

                %%%%% make save path if not already created %%%%%
                if ~isfolder(decResultsDir)
                    mkdir(decResultsDir)
                end

                %load data
                data = load(dayDatafname);
                data = data.data;
                if id == 2 || id == 3%need lap data for trial and SWR decoding
                    lapData = load(fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\JK' num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' params.decoding.decID{1} 'Data.mat']));%variable name is "lapData"
                    lapData = lapData.data;
                end%if id == 2 || id == 3

                %loop through environments
                for ee = 1:2%:length(params.environments)
                    currEnv = params.environments{ee};

                    %only loop if session exists (included because og sessions
                    % have fewer entries than update sessions)
                    if ~isempty(data.(currEnv))

                        %determine method to use
                        if strcmp(params.decoding.analysisType, 'kfold')
                            includedTrials = 0;
                            %find number of included trials
                            if strcmp(params.decoding.decID{id}, 'lap')
                                includedTrials = find(data.(currEnv).degBinOccup(end,:,1)~=0);
                            elseif strcmp(params.decoding.decID{id}, 'trial')
                                includedTrials = find(data.(currEnv).degBinOccup(end,end,:,1)~=0);
                            end%strcmp

                            %ensure enough trials
                            if length(includedTrials) >= params.decoding.minTrials
                                indices = crossvalind('Kfold', length(includedTrials), params.decoding.kfold);
                                for k = 1:params.decoding.kfold
                                    %%%%% get training/testing trials %%%%%
                                    trainInd = find(indices ~= k);
                                    testInd = find(indices == k);
                                    trainTrialInd = includedTrials(trainInd);
                                    testTrialInd = includedTrials(testInd);

                                    if strcmp(params.decoding.decID{id}, 'lap')
                                        %get training data
                                        trainingData = [];
                                        tmpDegEdges = data.(currEnv).degBinEdges;
                                        tmpDegEdges(1) = 0.0001;
                                        for clu = 1:length(data.(currEnv).apData(end).ID)
                                            %create ratemap: collect spikes per spatial bin
                                            % across all trials, sum occupancy per spatial bin
                                            % across trials, smooth each, then divide smoothed spike
                                            % counts by smoothed occupancy
                                            tmpSpikeCounts = []; tmpOccup = []; tmpRateMap = [];
                                            tmpSpikeCounts = gaussSmooth(histcounts(squeeze(data.(currEnv).spikePosAll(end,trainTrialInd,:,clu)), tmpDegEdges), 5);
                                            tmpOccup = gaussSmooth(squeeze(nansum(data.(currEnv).degBinOccup(end, trainTrialInd, :), 2))', 5);
                                            trainingData(clu,:) = tmpSpikeCounts ./ tmpOccup;%cell x position bin
                                        end%clu
                                        trainingData = trainingData.*(params.binsize_mstime/1000);%expected counts in the testing time bin size in s
                                        %change bins where cells don't fire to small value
                                        for clu = 1:size(trainingData,1)
                                            if sum(trainingData(clu,:) == 0) > 0
                                                zeroInd = trainingData(clu,:) == 0;
                                                trainingData(clu,zeroInd) = 1e-15;
                                            end%if
                                        end%clu
               
                                        %loop through test trials
                                        allTestTrEstPos = [];
                                        allTestTrActualPos = [];
                                        allTestTrActualPosDegBin = [];
                                        allTestTrError = [];
                                        for tr = 1:length(testTrialInd)
                                            %get testing data
                                            testingData = [];
                                            tmpTBinSpikeCount = [];
                                            tmpSpikeTime = squeeze(data.(currEnv).spikeTimeAll(end,testTrialInd(tr),:,1:length(data.(currEnv).apData(end).ID)));
                                            tmpTestSpikePos = squeeze(data.(currEnv).spikePosAll(end,testTrialInd(tr),:,1:length(data.(currEnv).apData(end).ID)));
                                            tmpActualPos = squeeze(data.(currEnv).tBinActualPos(end,testTrialInd(tr),:));
                                            tmpActualPos(tmpActualPos==0) = [];%cut down to size of this trial
                                            tmpTBins = squeeze(data.(currEnv).tBinEdges(end,testTrialInd(tr),1:length(tmpActualPos)));
                                            tmpSpeed = squeeze(data.(currEnv).tBinSpeed(end,testTrialInd(tr),1:length(tmpActualPos)));
                                            tmpTBins(tmpSpeed<params.speedTh) = [];%remove bins below moving threshold
                                            tmpActualPos(tmpSpeed<params.speedTh) = [];%remove bins below moving threshold
                                            for clu = 1:length(data.(currEnv).apData(end).ID)
                                                tmpTBinSpikeCount(clu,:) = histcounts(tmpSpikeTime(tmpSpikeTime(:,clu)~=0,clu),tmpTBins);
                                            end
                                            %output = cell x time bin
                                            testingData = tmpTBinSpikeCount;

                                            %%%%% DO DECODING %%%%%
                                            %adapted from getspatialprob.m
                                            spatialProb = [];
                                            totalSpikes = [];
                                            isLowOcc = isnan(trainingData(1,:));
                                            for b = 1:size(testingData,2)
                                                binSpikeCount = repmat(testingData(:,b), 1, size(trainingData,2));
                                                prob = prod(((trainingData.^binSpikeCount)./gamma(binSpikeCount + 1)).*exp(-trainingData),1,"omitmissing")'; % XZ modified 2/23/2025

                                                %if trainingCounts(:,b) is NaN, that means the occupancy was too low
                                                %in bin b. This should be the same across all cells. make the
                                                %spatial Prob 0 in those cases
                                                prob(isLowOcc,:) = 0;
                                                if sum(prob) ~= 0  % added by xz 02162025.
                                                    normProb = prob/sum(prob); % normalize across space to make the probabilities add up to 1
                                                else
                                                    normProb = prob;
                                                end
                                                if all(prob == 0)
                                                    normProb = ones(size(prob)) / length(prob); % Assign uniform probability instead of zeros
                                                end

                                                spatialProb(:,b) = normProb;
                                                totalSpikes(:,b) = sum(testingData(:,b));
                                            end%b

                                            %get decoding error for each time bin
                                            spatialProb = spatialProb(:,totalSpikes>2);%remove time bins with less than 2 spikes
                                            tmpActualPos = tmpActualPos(totalSpikes>2)';%remove time bins with less than 2 spikes
                                            [~, iMax] = max(spatialProb, [], 1);%find max position estimate per time bin
                                            tmpError = data.(currEnv).degBinEdges(iMax) - tmpActualPos;

                                            %save to output structs for all test trials
                                            tmpConcatInd = size(allTestTrEstPos,2);
                                            allTestTrError(:,tmpConcatInd+1:tmpConcatInd+size(tmpError,2)) = tmpError;
                                            allTestTrEstPos(:,tmpConcatInd+1:tmpConcatInd+size(spatialProb,2)) = spatialProb;
                                            allTestTrActualPos(:,tmpConcatInd+1:tmpConcatInd+size(tmpActualPos,2)) = tmpActualPos;
                                            allTestTrActualPosDegBin(:,tmpConcatInd+1:tmpConcatInd+size(tmpActualPos,2)) = discretize(tmpActualPos,tmpDegEdges);

                                        end%tr

                                        %get decoding estimate average across positions and trials
                                        allTestTrEstPosMn = [];
                                        allTestTrErrorMn = [];
                                        for b = 1:size(allTestTrEstPos,1)
                                            %ESTIMATED POSITION%
                                            % output = estimated position x acutal position
                                            allTestTrEstPosMn(:,b) = nanmean(allTestTrEstPos(:,allTestTrActualPosDegBin==b),2);
                                            %ERROR%
                                            allTestTrErrorMn(b) = nanmean(allTestTrError(:,allTestTrActualPosDegBin==b),2);
                                        end%b

                                        %add to output structs
                                        decodedEstPosAll(ee,ss,k,:,:) = allTestTrEstPosMn;
                                        decodedErrorAll(ee,ss,k,:) = allTestTrErrorMn;

                                    elseif strcmp(params.decoding.decID{id}, 'trial')
                                        trainingData = [];
                                        tmpDegBinSpikeCount = [];
                                        tmpDegBinSpikeCountSmooth = [];
                                        tmpDegBinRateMap = [];
                                        for znType = 1:size(data.(currEnv).degBinEdges,2)
                                            for t = 1:length(trainTrialInd)
                                                tmpSpikePos = squeeze(data.(currEnv).spikePosAll(ss,znType,trainTrialInd(t),:,1:length(data.(currEnv).apData(ss).ID)));
                                                tmpOccupSmooth = transpose(squeeze(data.(currEnv).degBinOccupSmooth(ss,znType,trainTrialInd(t),:)));
                                                %in case trial wraps around 360, put pos bins into ascending order for histcounts below
                                                [tmpTrainDegBins, tmpTrainDegBinsInd] = sort(squeeze(data.(currEnv).degBinEdges(ss,znType,trainTrialInd(t),:)));
                                                for clu = 1:length(data.(currEnv).apData(ss).ID)
                                                    tmpDegBinSpikeCount(clu,:) = histcounts(tmpSpikePos(tmpSpikePos(:,clu)~=0,clu),tmpTrainDegBins');
                                                    [~, tmpResortedTrainDegBinsInd] = sort(tmpTrainDegBinsInd(tmpTrainDegBinsInd~=length(tmpTrainDegBinsInd)));%all bins except last
                                                    tmpDegBinSpikeCount(clu,:) = tmpDegBinSpikeCount(clu,tmpResortedTrainDegBinsInd);%return to original order
                                                    tmpDegBinSpikeCountSmooth(clu,:) = gaussSmooth(tmpDegBinSpikeCount(clu,:),2);
                                                    tmpDegBinRateMap(clu,:) = tmpDegBinSpikeCountSmooth(clu,:) ./ tmpOccupSmooth;
                                                end
                                                %output = trial x zntype x cell x position bin
                                                trainingData(t,znType,:,:) = tmpDegBinRateMap;
                                            end%znType
                                        end%t
                                        trainingData = squeeze(nanmean(trainingData,1));%average across laps: znType x cell x position
                                        trainingData = trainingData.*(params.binsize_mstime/1000);%expected counts in the testing time bin size in s

                                        %loop through trial types
                                        for znType = 1:size(data.(currEnv).degBinEdges,2)%1=reward, 2=nonreward, 3 = alt nonreward
                                            %get training data for this zone type
                                            tmpTrainingData = [];
                                            tmpTrainingData = squeeze(trainingData(znType,:,:));%cell x position
                                            %change bins where cells don't fire to small value
                                            for clu = 1:size(tmpTrainingData,1)
                                                if sum(tmpTrainingData(clu,:) == 0) > 0
                                                    zeroInd = tmpTrainingData(clu,:) == 0;
                                                    tmpTrainingData(clu,zeroInd) = 1e-15;
                                                end%if
                                            end%clu

                                            %loop through test trials
                                            allTestTrEstPos = [];
                                            allTestTrActualPos = [];
                                            allTestTrActualPosDegBin = [];
                                            allTestTrError = [];
                                            for tr = 1:length(testTrialInd)
                                                %get testing data
                                                testingData = [];
                                                tmpTBinSpikeCount = [];
                                                tmpSpikeTime = squeeze(data.(currEnv).spikeTimeAll(ss,znType,testTrialInd(tr),:,1:length(data.(currEnv).apData(ss).ID)));
                                                tmpTestSpikePos = squeeze(data.(currEnv).spikePosAll(ss,znType,testTrialInd(tr),:,1:length(data.(currEnv).apData(ss).ID)));
                                                tmpActualPos = squeeze(data.(currEnv).tBinActualPos(ss,znType,testTrialInd(tr),:));
                                                tmpActualPos(tmpActualPos==0) = [];%cut down to size of this trial
                                                tmpTBins = squeeze(data.(currEnv).tBinEdges(ss,znType,testTrialInd(tr),1:length(tmpActualPos)));
                                                tmpSpeed = squeeze(data.(currEnv).tBinSpeed(ss,znType,testTrialInd(tr),1:length(tmpActualPos)));
                                                tmpTBins(tmpSpeed<params.speedTh) = [];%remove bins below moving threshold
                                                tmpActualPos(tmpSpeed<params.speedTh) = [];%remove bins below moving threshold
                                                for clu = 1:length(data.(currEnv).apData(ss).ID)
                                                    tmpTBinSpikeCount(clu,:) = histcounts(tmpSpikeTime(tmpSpikeTime(:,clu)~=0,clu),tmpTBins);
                                                end
                                                %output = cell x time bin
                                                testingData = tmpTBinSpikeCount;

                                                %%%%% DO DECODING %%%%%
                                                %adapted from getspatialprob.m
                                                spatialProb = [];
                                                totalSpikes = [];
                                                isLowOcc = isnan(tmpTrainingData(1,:));
                                                for b = 1:size(testingData,2)
                                                    binSpikeCount = repmat(testingData(:,b), 1, size(tmpTrainingData,2));
                                                    prob = prod(((tmpTrainingData.^binSpikeCount)./gamma(binSpikeCount + 1)).*exp(-tmpTrainingData),1,"omitmissing")'; % XZ modified 2/23/2025

                                                    %if trainingCounts(:,b) is NaN, that means the occupancy was too low
                                                    %in bin b. This should be the same across all cells. make the
                                                    %spatial Prob 0 in those cases
                                                    prob(isLowOcc,:) = 0;
                                                    if sum(prob) ~= 0  % added by xz 02162025.
                                                        normProb = prob/sum(prob); % normalize across space to make the probabilities add up to 1
                                                    else
                                                        normProb = prob;
                                                    end
                                                    if all(prob == 0)
                                                        normProb = ones(size(prob)) / length(prob); % Assign uniform probability instead of zeros
                                                    end

                                                    spatialProb(:,b) = normProb;
                                                    totalSpikes(:,b) = sum(testingData(:,b));
                                                end%b

                                                %get decoding error for each time bin
                                                spatialProb = spatialProb(:,totalSpikes>2);%remove time bins with less than 2 spikes
                                                tmpActualPos = tmpActualPos(totalSpikes>2)';%remove time bins with less than 2 spikes
                                                [~, iMax] = max(spatialProb, [], 1);%find max position estimate per time bin
                                                tmpUnsortedTestDegBins = squeeze(data.(currEnv).degBinEdges(ss,znType,testTrialInd(tr),1:size(data.(currEnv).degBinEdges,4)-1));
                                                tmpError = tmpUnsortedTestDegBins(iMax)' - tmpActualPos;
                                                tmpError(tmpError > 180) = tmpError(tmpError > 180) - 360;
                                                tmpError(tmpError < -180) = 360 + tmpError(tmpError < -180);

                                                %save to output structs for all test trials
                                                tmpConcatInd = size(allTestTrEstPos,2);
                                                allTestTrError(:,tmpConcatInd+1:tmpConcatInd+size(tmpError,2)) = tmpError;
                                                allTestTrEstPos(:,tmpConcatInd+1:tmpConcatInd+size(spatialProb,2)) = spatialProb;
                                                allTestTrActualPos(:,tmpConcatInd+1:tmpConcatInd+size(tmpActualPos,2)) = tmpActualPos;
                                                %in case trial wraps around 360, put pos bins into ascending order for discretize below
                                                [tmpTestDegBins, tmpTestDegBinsInd] = sort(squeeze(data.(currEnv).degBinEdges(ss,znType,testTrialInd(tr),:))');
                                                allTestTrActualPosDegBin(:,tmpConcatInd+1:tmpConcatInd+size(tmpActualPos,2)) = discretize(tmpActualPos,tmpTestDegBins);
                                                [~, tmpResortedTestDegBinsInd] = sort(tmpTrainDegBinsInd(tmpTrainDegBinsInd~=length(tmpTrainDegBinsInd)));%all bins except last
                                                tmpTestDegBinsInd = tmpTestDegBinsInd(tmpTestDegBinsInd~=length(tmpTestDegBinsInd));%%all bins except last
                                                for tmpInd = 1:length(tmpResortedTestDegBinsInd)
                                                    allTestTrActualPosDegBin(allTestTrActualPosDegBin == tmpResortedTestDegBinsInd(tmpInd)) = tmpTestDegBinsInd(tmpInd);
                                                end%tmpInd
                                            end%tr

                                            %get decoding estimate average across positions and trials
                                            allTestTrEstPosMn = [];
                                            allTestTrErrorMn = [];
                                            for b = 1:size(allTestTrEstPos,1)
                                                %ESTIMATED POSITION%
                                                % output = estimated position x acutal position
                                                allTestTrEstPosMn(:,b) = nanmean(allTestTrEstPos(:,allTestTrActualPosDegBin==b),2);
                                                %ERROR%
                                                allTestTrErrorMn(b) = nanmean(allTestTrError(:,allTestTrActualPosDegBin==b),2);
                                            end%b

                                            %add to output structs
                                            decodedEstPosAll(ee,ss,k,znType,:,:) = allTestTrEstPosMn;
                                            decodedErrorAll(ee,ss,k,znType,:) = allTestTrErrorMn;
                                        end%znType
                                    end%trials or laps
                                end%k

                                %plot and save fig
                                if strcmp(params.decoding.decID{id}, 'lap')%for laps
                                    %plot
                                    figure('Position',[0 0 1500 1500]);
                                    %clims = [0 0.75*max(decOutMatEstPos(ee,ss,:,:,:),[],'all')];
                                    %imagesc(1:size(decOutMatEstPos,5), 1:size(decOutMatEstPos,5), squeeze(nanmean(decOutMatEstPos(ee,ss,:,:,:),3)), clims);
                                    imagesc(1:size(decodedEstPosAll,5), 1:size(decodedEstPosAll,5), squeeze(nanmean(decodedEstPosAll(ee,ss,:,:,:),3)));
                                    colormap('hot'); colorbar;
                                    title(sprintf('%s%d_%d_%s', params.iden, data.(currEnv).sessionInfo{end}(1,1), data.(currEnv).sessionInfo{end}(1,2), currEnv), Interpreter="none")
                                    %save
                                    figdir = [dirs.savefigures '\Neural\Decoding\Lap\' sprintf('JK%d_%d', data.(currEnv).sessionInfo{end}(1,1), data.(currEnv).sessionInfo{end}(1,2))];
                                    if ~isfolder(figdir)
                                        mkdir(figdir)
                                    end
                                    filename = [figdir '\' sprintf('decoding_kfold_laps_tbins')];
                                    print(gcf,filename,'-dpng','-r300')
                                elseif strcmp(params.decoding.decID{id}, 'trial')%for trials
                                    for znType = 1:size(trainingData,1)
                                        if znType == 1
                                            znTypeNm = 'RZ';
                                        elseif znType == 2
                                            znTypeNm = 'NRZ';
                                        elseif znType == 3
                                            znTypeNm = 'NevRZ';
                                        end%znTypeNm
                                        %plot
                                        figure('Position',[0 0 1500 1500]);
                                        imagesc(1:size(decodedEstPosAll,6), 1:size(decodedEstPosAll,6), squeeze(nanmean(decodedEstPosAll(ee,ss,:,znType,:,:),3)));
                                        colormap('hot'); colorbar;
                                        title(sprintf('%s%d_%d_%s_%s', params.iden, data.(currEnv).sessionInfo{ss}(1,1), data.(currEnv).sessionInfo{ss}(1,2), currEnv, znTypeNm), Interpreter="none")
                                        %save
                                        figdir = [dirs.savefigures '\Neural\Decoding\Trial\' sprintf('JK%d_%d', data.(currEnv).sessionInfo{ss}(1,1), data.(currEnv).sessionInfo{ss}(1,2))];
                                        if ~isfolder(figdir)
                                            mkdir(figdir)
                                        end
                                        filename = [figdir '\' sprintf('decoding_kfold_%s_trials_tbins', znTypeNm)];
                                        print(gcf,filename,'-dpng','-r300')
                                    end%znType
                                end%strcmp

                            elseif length(includedTrials) < params.decoding.minTrials && ~isempty(data.(currEnv).sessionInfo{ss})
                                sprintf('Not enough included %s %ss (<%d) for JK%d_%d, only %d.', currEnv, params.decoding.decID{id}, ...
                                    params.decoding.minTrials, data.(currEnv).sessionInfo{ss}(1,1), data.(currEnv).sessionInfo{ss}(1,2), length(includedTrials))
                            else%did not do this environment
                            end%included trials > minTrials

                        elseif strcmp(params.decoding.analysisType, 'leave1out')
                            %Note: Currently built for only trial and ripple analyses

                            if strcmp(params.decoding.decID{id}, 'trial')
                                %loop through trial types
                                for znType = 1:size(data.(currEnv).degBinEdges,2)%1=reward, 2=nonreward, 3 = alt nonreward
                                    %loop through trials, cutting down to only completed laps
                                    lpCtr = 0;
                                    numCmpltLap = sum(sum(lapData.(currEnv).degBinOccupSmooth(ss,:,:),3)~=0);
                                    numCmpltLapTr = find(data.(currEnv).degBinOccupSmooth(ss,end,:,1)~=0,1,'last');
                                    for tr = 1:numCmpltLapTr
                                        %add one to counter each start of lap
                                        if any(tr == 1:length(data.(currEnv).azBins_deg):numCmpltLapTr)
                                            lpCtr = lpCtr + 1;
                                        end%lpCtr

                                        %get testing data
                                        testingData = [];
                                        tmpTBinSpikeCount = [];
                                        tmpTestSpikeTime = squeeze(data.(currEnv).spikeTimeAll(ss,znType,tr,:,1:length(data.(currEnv).apData(ss).ID)));
                                        tmpActualPos = squeeze(data.(currEnv).tBinActualPos(ss,znType,tr,:));
                                        tmpActualPos(tmpActualPos==0) = [];%cut down to size of this trial
                                        tmpTBins = squeeze(data.(currEnv).tBinEdges(ss,znType,tr,1:length(tmpActualPos)));
                                        tmpSpeed = squeeze(data.(currEnv).tBinSpeed(ss,znType,tr,1:length(tmpActualPos)));
                                        tmpTBins(tmpSpeed<params.speedTh) = [];%remove bins below moving threshold
                                        tmpActualPos(tmpSpeed<params.speedTh) = [];%remove bins below moving threshold
                                        for clu = 1:length(data.(currEnv).apData(ss).ID)
                                            tmpTBinSpikeCount(clu,:) = histcounts(tmpTestSpikeTime(tmpTestSpikeTime(:,clu)~=0,clu),tmpTBins);
                                        end
                                        %output = cell x time bin
                                        testingData = tmpTBinSpikeCount;

                                        %get training data
                                        trainingData = [];
                                        tmpDegBinSpikeCount = [];
                                        tmpDegBinSpikeCountSmooth = [];
                                        tmpDegBinRateMap = [];
                                        tmpTrainDegBins = lapData.(currEnv).degBinEdges;
                                        for t = 1:numCmpltLap
                                            tmpOccupSmooth = transpose(squeeze(lapData.(currEnv).degBinOccupSmooth(ss,t,:)));
                                            for clu = 1:length(lapData.(currEnv).apData(ss).ID)
                                                tmpTrainSpikeTime = squeeze(lapData.(currEnv).spikeTimeAll(ss,t,:,clu));
                                                tmpTrainSpikePos = squeeze(lapData.(currEnv).spikePosAll(ss,t,:,clu));
                                                tmpTrainSpikePos(find(ismember(tmpTrainSpikeTime, tmpTestSpikeTime(tmpTestSpikeTime(:,clu)~=0,clu)))) = [];%remove test spikes from training data
                                                tmpDegBinSpikeCount(clu,:) = histcounts(tmpTrainSpikePos(tmpTrainSpikePos~=0),tmpTrainDegBins');
                                                tmpDegBinSpikeCountSmooth(clu,:) = gaussSmooth(tmpDegBinSpikeCount(clu,:),2);
                                                tmpDegBinRateMap(clu,:) = tmpDegBinSpikeCountSmooth(clu,:) ./ tmpOccupSmooth;
                                            end
                                            %output = laps x cell x position bin
                                            trainingData(t,:,:) = tmpDegBinRateMap;
                                        end%t
                                        trainingData = squeeze(nanmean(trainingData,1));%average across laps: cell x position
                                        trainingData = trainingData.*(params.binsize_mstime/1000);%expected counts in the testing time bin size in s
                                        %change bins where cells don't fire to small value
                                        for clu = 1:size(trainingData,1)
                                            if sum(trainingData(clu,:) == 0) > 0
                                                zeroInd = trainingData(clu,:) == 0;
                                                trainingData(clu,zeroInd) = 1e-15;
                                            end%if
                                        end%clu

                                        %%%%% DO DECODING %%%%%
                                        %adapted from getspatialprob.m
                                        spatialProb = [];
                                        totalSpikes = [];
                                        isLowOcc = isnan(trainingData(1,:));
                                        for b = 1:size(testingData,2)
                                            binSpikeCount = repmat(testingData(:,b), 1, size(trainingData,2));
                                            prob = prod(((trainingData.^binSpikeCount)./gamma(binSpikeCount + 1)).*exp(-trainingData),1,"omitmissing")'; % XZ modified 2/23/2025

                                            %if trainingCounts(:,b) is NaN, that means the occupancy was too low
                                            %in bin b. This should be the same across all cells. make the
                                            %spatial Prob 0 in those cases
                                            prob(isLowOcc,:) = 0;
                                            if sum(prob) ~= 0  % added by xz 02162025.
                                                normProb = prob/sum(prob); % normalize across space to make the probabilities add up to 1
                                            else
                                                normProb = prob;
                                            end
                                            if all(prob == 0)
                                                normProb = ones(size(prob)) / length(prob); % Assign uniform probability instead of zeros
                                            end

                                            spatialProb(:,b) = normProb;
                                            totalSpikes(:,b) = sum(testingData(:,b));
                                        end%b

                                        %get decoding error for each time bin
                                        spatialProb = spatialProb(:,totalSpikes>2);%remove time bins with less than 2 spikes
                                        tmpActualPos = tmpActualPos(totalSpikes>2)';%remove time bins with less than 2 spikes
                                        [~, iMax] = max(spatialProb, [], 1);%find max position estimate per time bin
                                        tmpError = tmpTrainDegBins(iMax) - tmpActualPos;
                                        tmpError(tmpError > 180) = tmpError(tmpError > 180) - 360;
                                        tmpError(tmpError < -180) = 360 + tmpError(tmpError < -180);

                                        %get decoding estimate average across positions
                                        tmpTestTrActualPosDegBin = [];
                                        tmpTestTrEstPosMn = [];
                                        tmpTestTrErrorMn = [];
                                        tmpTestTrActualPosDegBin = discretize(tmpActualPos,tmpTrainDegBins);

                                        for b = 1:size(spatialProb,1)
                                            %ESTIMATED POSITION%
                                            % output = estimated position x acutal position
                                            tmpTestTrEstPosMn(:,b) = nanmean(spatialProb(:,tmpTestTrActualPosDegBin==b),2);
                                            %ERROR%
                                            tmpTestTrErrorMn(b) = nanmean(tmpError(:,tmpTestTrActualPosDegBin==b),2);
                                        end%b

                                        %save to output structs
                                        decodedEstPosAll(ee,ss,znType,tr,:,:) = tmpTestTrEstPosMn;
                                        decodedErrorAll(ee,ss,znType,tr,:) = tmpTestTrErrorMn;
                                    end%tr
                                end%znType

                                % % % %save data for this session
                                % % % save(decResultsfname

                            elseif strcmp(params.decoding.decID{id}, 'ripple')
                                %Note: currently 1 position bin per ripple (JLK 1/19/26)
                                for ch = 1:size(data.(currEnv).degBinRateMap,4)
                                    for r = 1:size(data.(currEnv).degBinRateMap,5)
                                        %get testing data
                                        testingData = squeeze(data.(currEnv).degBinRateMap(ss,znType,tr,ch,r,:));%cell x position
                                        % % %normalize per unit
                                        % % spikingNeuronidx = min(testingData,[],2) ~= max(testingData,[],2);
                                        % % testingData(spikingNeuronidx, :) = (testingData(spikingNeuronidx, :) - min(testingData(spikingNeuronidx, :))) ./ (max(testingData(spikingNeuronidx, :)) - min(testingData(spikingNeuronidx, :)));
                                        % % testingData(~spikingNeuronidx, :) = zeros(sum(~spikingNeuronidx), size(testingData, 2));
                                        testingData = testingData .* .001;%params.decoding.scaleFactor;%scale for plotting
                                        % % %change bins where cells don't fire to small value
                                        % % for clu = 1:size(testingData,1)
                                        % %     if sum(testingData(clu,:) == 0) > 0
                                        % %         zeroInd = testingData(clu,:) == 0;
                                        % %         testingData(clu,zeroInd) = 1e-15;
                                        % %     end%if
                                        % % end%clu
                                        if size(testingData,1)==1; testingData = testingData';end

                                        %%%%% DO DECODING %%%%%
                                        %adapted from getspatialprob.m
                                        spatialProb = [];
                                        isLowOcc = isnan(tmpTrainingData(1,:));
                                        for b = 1:size(testingData,2)
                                            binSpikeCount = repmat(testingData(:,b), 1, size(tmpTrainingData,2));
                                            prob = prod(((tmpTrainingData.^binSpikeCount)./gamma(binSpikeCount + 1)).*exp(-tmpTrainingData),1,"omitmissing")'; % XZ modified 2/23/2025

                                            %if trainingCounts(:,b) is NaN, that means the occupancy was too low
                                            %in bin b. This should be the same across all cells. make the
                                            %spatial Prob 0 in those cases
                                            prob(isLowOcc,:) = 0;
                                            if sum(prob) ~= 0  % added by xz 02162025.
                                                normProb = prob/sum(prob); % normalize across space to make the probabilities add up to 1
                                            else
                                                normProb = prob;
                                            end
                                            if all(prob == 0)
                                                normProb = ones(size(prob)) / length(prob); % Assign uniform probability instead of zeros
                                            end

                                            spatialProb(:,b) = normProb;
                                        end%b

                                        %save to output structs
                                        decOut(ee,ss,znType,tr,ch,r,:).spatialProb = spatialProb; %bins x time
                                        decOutMat(ee,ss,znType,tr,ch,r,:,:) = spatialProb;

                                    end%r
                                end%ch
                            end%strcmp
                            %     end%lp
                            % end%znType
                        end%params.decoding.analysisType

                    end%if data exists
                end%ee           
        end%if session not done
    end%ss

    %save data
    decodedEstPosAllfname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\' params.decoding.decID{id} 'DecodedEstPosAll.mat']);
    save(decodedEstPosAllfname, 'decodedEstPosAll', '-v7.3');
    decodedErrorAllfname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\' params.decoding.decID{id} 'DecodedErrorAll.mat']);
    save(decodedErrorAllfname, 'decodedErrorAll', '-v7.3');

end%id

    %%%%%%%%%%%%%%%%
    %%%%% PLOT %%%%%
    %%%%%%%%%%%%%%%%

    %%%%% Trials %%%%%
    for ee = 1:2
        currEnv = params.environments{ee};
        for ss = 1:3%:size(data.(currEnv).degBinEdges,1)
            if size(data.(currEnv).degBinOccupSmooth,2)>=ss%if session exists
                for znType = 1:size(data.(currEnv).degBinEdges,2)
                    %INDIVIDUAL TRIALS%
                    numCmpltLapTr = find(data.(currEnv).degBinOccupSmooth(ss,end,:,1)~=0,1,'last');
                    for tr = 1:numCmpltLapTr
                        %plot decoding results
                        figure; hold on
                        tmpData = squeeze(decodedEstPosAll(ee,ss,znType,tr,:,:));
                        tmpDataXVals = find(~isnan(tmpData(1,:)),1,'first'):find(~isnan(tmpData(1,:)),1,'last');
                        tmpCurrZn = squeeze(data.(currEnv).degBinEdges(ss,znType,tr,1:size(data.(currEnv).degBinEdges,4)))/2;
                        if tmpDataXVals(end)~= tmpCurrZn(end)
                            tmpDataXVals(end+1) = tmpCurrZn(end);
                        elseif tmpDataXVals(1)~= tmpCurrZn(2)
                            tmpDataXVals = [NaN tmpDataXVals];
                        end%fix length of tmpDataXVals
                        tmpPlotXVals = tmpDataXVals(1)-.5:tmpDataXVals(end)+.5;
                        imagesc(tmpDataXVals, 1:size(tmpData,1), tmpData(:,tmpDataXVals));
                        colormap('hot'); colorbar;
                        xlim([tmpPlotXVals(1) tmpPlotXVals(end)])
                        ylim([1 size(tmpData,1)])
                        if ee == 1%original
                            plotColors = [1 0 0; .9 .9 .9];
                        elseif ee == 2%update
                            plotColors = [0 1 0; 1 0 0; .9 .9 .9];
                        end

                        %draw lines for zones
                        %DECODED REWARD ZONES%
                        for znNum = 1:length(data.(currEnv).azBins_deg)
                            %RZ
                            plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).azBins_deg(znNum)+10)/2,'--', 'Color', plotColors(1,:))
                            plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).azBins_deg(znNum)+20)/2,'--', 'Color', plotColors(1,:))
                        end%znNum RZ
                        %DECODED NONREWARD ZONES%
                        for znNum = 1:length(data.(currEnv).czBins_deg)
                            %NRZ
                            plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).czBins_deg(znNum))/2,'--', 'Color', plotColors(2,:))
                            plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).czBins_deg(znNum)+10)/2,'--', 'Color', plotColors(2,:))
                        end%znNum NRZ
                        if ee == 2%Update
                            %DECODED NEVER REWARDED ZONES%
                            for znNum = 1:length(data.(currEnv).nevrzBins_deg)
                                %NevRZ
                                plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).nevrzBins_deg(znNum))/2,'--', 'Color', plotColors(3,:))
                                plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).nevrzBins_deg(znNum)+10)/2,'--', 'Color', plotColors(3,:))
                            end%znNum NevRZ
                        end%ee
                        %DECODED CURRENT ZONE%
                        if ee == 1%original track
                            sessionType = 'original';
                            vline(tmpPlotXVals(11),'r--')
                            vline(tmpPlotXVals(16),'r--')
                        elseif ee == 2%update track
                            sessionType = 'update';
                            vline(tmpPlotXVals(11),'g--')
                            vline(tmpPlotXVals(16),'g--')
                        end%ee
                        if any(abs(diff(tmpCurrZn))>100) == 0%no wrap arround 360 deg
                            plot(tmpPlotXVals,tmpCurrZn,'--y')
                        elseif any(abs(diff(tmpCurrZn))>100) == 1%wrap arround 360 deg
                            tmpCurrZnLine1 = tmpCurrZn;
                            tmpCurrZnLine1(tmpCurrZnLine1<100) = nan;
                            plot(tmpPlotXVals,tmpCurrZnLine1,'--y')
                            tmpCurrZnLine2 = tmpCurrZn;
                            tmpCurrZnLine2(tmpCurrZnLine2>100) = nan;
                            plot(tmpPlotXVals,tmpCurrZnLine2,'--y')
                        end%any

                        %save
                        figdir = [dirs.savefigures '\Neural\Decoding\Trial\' sprintf('JK%d_%d', data.(currEnv).sessionInfo{ss}(1,1), data.(currEnv).sessionInfo{ss}(1,2))];
                        if ~isfolder(figdir)
                            mkdir(figdir)
                        end
                        if znType == 1
                            znTypeNm = 'RZ';
                        elseif znType == 2
                            znTypeNm = 'NRZ';
                        elseif znType == 3
                            znTypeNm = 'NevRZ';
                        end%znTypeNm
                        filename = [figdir '\' sprintf('decoding_%s_individual_%s_trials_%d_tbins', sessionType, znTypeNm, tr)];
                        print(gcf,filename,'-dpng','-r300')
                    end%tr

                    %AVERAGE ACROSS GIVEN ZONE%
                    for znNumTr = (1:length(data.(currEnv).azBins_deg))-1
                        %plot decoding results
                        tr = 1+znNumTr:length(data.(currEnv).azBins_deg):numCmpltLapTr;
                        if ~isempty(tr)
                            figure('Position',[0 0 1500 1500]); hold on
                            tmpData = squeeze(nanmean(decodedEstPosAll(ee,ss,znType,tr,:,:),4));
                            tmpDataXVals = find(~isnan(tmpData(1,:)),1,'first'):find(~isnan(tmpData(1,:)),1,'last');
                            tmpPlotXVals = tmpDataXVals(1)-.5:tmpDataXVals(end)+.5;
                            imagesc(tmpDataXVals, 1:size(tmpData,1), tmpData(:,tmpDataXVals));
                            colormap('hot'); colorbar;
                            xlim([tmpPlotXVals(1) tmpPlotXVals(end)])
                            ylim([1 size(tmpData,1)])
                            if ee == 1%original
                                plotColors = [1 0 0; .9 .9 .9];
                            elseif ee == 2%update
                                plotColors = [0 1 0; 1 0 0; .9 .9 .9];
                            end

                            %draw lines for zones
                            %DECODED REWARD ZONES%
                            for znNum = 1:length(data.(currEnv).azBins_deg)
                                %RZ
                                plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).azBins_deg(znNum)+10)/2,'--', 'Color', plotColors(1,:))
                                plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).azBins_deg(znNum)+20)/2,'--', 'Color', plotColors(1,:))
                            end%znNum RZ
                            %DECODED NONREWARD ZONES%
                            for znNum = 1:length(data.(currEnv).czBins_deg)
                                %NRZ
                                plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).czBins_deg(znNum))/2,'--', 'Color', plotColors(2,:))
                                plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).czBins_deg(znNum)+10)/2,'--', 'Color', plotColors(2,:))
                            end%znNum NRZ
                            if ee == 2%Update
                                %DECODED NEVER REWARDED ZONES%
                                for znNum = 1:length(data.(currEnv).nevrzBins_deg)
                                    %NevRZ
                                    plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).nevrzBins_deg(znNum))/2,'--', 'Color', plotColors(3,:))
                                    plot(tmpPlotXVals,ones(1,length(tmpPlotXVals))*(data.(currEnv).nevrzBins_deg(znNum)+10)/2,'--', 'Color', plotColors(3,:))
                                end%znNum NevRZ
                            end%ee
                            %DECODED CURRENT ZONE%
                            if ee == 1%original track
                                sessionType = 'original';
                                vline(tmpPlotXVals(11),'r--')
                                vline(tmpPlotXVals(16),'r--')
                            elseif ee == 2%update track
                                sessionType = 'update';
                                vline(tmpPlotXVals(11),'g--')
                                vline(tmpPlotXVals(16),'g--')
                            end%ee
                            tmpCurrZn = squeeze(data.(currEnv).degBinEdges(ss,znType,tr,1:size(data.(currEnv).degBinEdges,4)))/2;
                            if any(abs(diff(tmpCurrZn))>100) == 0%no wrap arround 360 deg
                                plot(tmpPlotXVals,tmpCurrZn(1):tmpCurrZn(end),'--y')
                            elseif any(abs(diff(tmpCurrZn))>100) == 1%wrap arround 360 deg
                                tmpCurrZnLine1 = tmpCurrZn;
                                tmpCurrZnLine1(tmpCurrZnLine1<100) = nan;
                                plot(tmpPlotXVals,tmpCurrZnLine1,'--y')
                                tmpCurrZnLine2 = tmpCurrZn;
                                tmpCurrZnLine2(tmpCurrZnLine2>100) = nan;
                                plot(tmpPlotXVals,tmpCurrZnLine2,'--y')
                            end%any

                            %save
                            figdir = [dirs.savefigures '\Neural\Decoding\Trial\' sprintf('JK%d_%d', data.(currEnv).sessionInfo{ss}(1,1), data.(currEnv).sessionInfo{ss}(1,2))];
                            if ~isfolder(figdir)
                                mkdir(figdir)
                            end
                            if znType == 1
                                znTypeNm = 'RZ';
                            elseif znType == 2
                                znTypeNm = 'NRZ';
                            elseif znType == 3
                                znTypeNm = 'NevRZ';
                            end
                            filename = [figdir '\' sprintf('decoding_%s_average_%s_trials_%d_tbins', sessionType, znTypeNm, znNumTr+1)];
                            print(gcf,filename,'-dpng','-r300')
                        end%if ~isempty(tr)
                    end%znNumTr
                end%znType
            end%if session exists
        end%ss
    end%ee

    %%%%% Ripples %%%%%
    %plot individual ripples
    figure('Position',[0 0 1500 1500]); hold on
    imagesc(squeeze( decOutMat(1,1,1,6,2,1,:)));
    colormap('hot'); colorbar;
    %plot(0.5:1.5,ones(1,2)*ripLoc(rip),'--r')




    % % %save data
    % % dir2save = fullfile(dirs.saveoutputstructs, 'Data\Behavior\ROC', [params.iden num2str(animal)], ...
    % %     num2str(sessionInfo(1,2)));
    % % if ~isfolder(dir2save); mkdir(dir2save); end
    % % save([dir2save, '\', params.rocID{id}, '.mat'], 'data', 'data');


else
    load(decfname);
end







end%function