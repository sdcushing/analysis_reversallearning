%%%%% add this section to getDirectories AndParams_JLK.m %%%%%
%decoding
params.decoding.decID = {'lap', 'trial'};
params.decoding.trialType = 'all'; %options: 'all', 'correct', 'incorrect'
params.decoding.analysisType = 'kfold'; %options: 'kfold', 'leave1out'
params.decoding.trainPer = 80;
params.decoding.testPer = 100-params.decoding.trainPer;
params.decoding.kfold = 100/params.decoding.testPer;
params.decoding.minTrials = params.decoding.kfold;
params.decoding.cellType = 'all'; %options: 'all', 'pyr', 'narrow'
params.decoding.cellLocation = 'all'; %options: 'all', 'CA1', 'CA3', 'HIP', 'cortical'
params.decoding.correctCriterion = 'reward'; %options: 'reward', 'anticipation'
params.decoding.scaleFactor = 4;%scaling factor to make values intergers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params.nDeg = 100; %360 for lap trials; 100 for reward center trials (30 deg before and 60 deg after RZ)
params.kfold = 5; %two folds cross validation, using 4/5 of data for training and the remaining 1/5 for testing
params.pfBinsize = 2; %in degrees, for making rate maps over space
params.decodeBinsize = 0.25; %decoding time binsize in sec
params.placefields = 0;
params.speedThreshold = 1; %in deg/s
params.posStep = 0.02; %in sec
params.occThreshold = 0.1; %in sec
params.cellCountThreshold = 1;
params.minCorrectTrials = 20;
params.minNCells = 20;
params.environments = {'Fam','Nov'};
params.zones = [1 2];
params.similarity = 0;
params.nonThetaRipples = 0; %1 for nontheta only ripples, 0 for all ripples
params.timeAroundMid = 0.125; %in s
params.decodingBin_s = params.timeAroundMid * 2; %use one timebin per ripple (i.e. do not break into multiple timebins)
% params.decodingBin_s = 0.020; %20ms decoding bins
params.decodingBin_samp = params.decodingBin_s .* params.ap_samprate;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [decOut] = getDecodingResults_linearJLK(allindex,dirs,uniqSess,params)
%adapted and integrated from position_decoding_function.m,
% runbayesiandecoding.m, decodecurrentposition.m/decodecurrentposition_alltrials.m,
% gettrainingdata.m/gettestingdata.m, and getpositionestimate.m
%code follows similar structure as getBehaviorROC_JLK.m, but only save data
% across all sessions, not per session and all sessions

%% create or load decoding data for each session %%
sessionInfo = [];

for id = 1:length(params.decoding.decID)

    decDatafname = fullfile(dirs.saveoutputstructs, 'Data\Neural\Decoding', [params.decoding.decID{id} 'Data.mat']);

    if ~exist(decDatafname) || params.rewrite.decodingData

        %%%%% make save path if not already created %%%%%
        if ~isfolder([dirs.saveoutputstructs, 'Data\Neural\Decoding'])
            mkdir([dirs.saveoutputstructs, 'Data\Neural\Decoding'])
        end

        %%%%% initialize data structure with empty fields %%%%%
        data = [];
        for ee = 1:length(params.environments)
            data.(params.environments{ee}) = [];
        end

        %loop through unique sessions
        for ss = 1:8%size(uniqSess,1)

            animal = uniqSess(ss, 1);
            sessDate = uniqSess(ss, 2);
            tempidx = find(ismember(allindex(:,1), animal) & ismember(allindex(:,2), sessDate) & ismember(allindex(:,6), 2) & ismember(allindex(:,5), 1:4));%include this animal, this date, active and recording only
            sessionInfo = allindex(tempidx, :);

            if ~isempty(sessionInfo)

                %loop through files to COLLECT DATA ACROSS SESSIONS FOR THIS DAY
                lpCtr = zeros(3,1);
                files = sessionInfo(:,3);
                for f = 1:length(files)
                    currFile = files(f);
                    sessionType = sessionInfo(f,4);

                    %%%%% load data %%%%%
                    %lap or trial data for decoding
                    if strcmp(params.decoding.decID{id}, 'lap')
                        neuralDataFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                            [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'rawDataByLapNeural.mat');
                        if isfile(neuralDataFname)
                            rawDataByLapNeural = load(neuralDataFname);
                            rawDataByLapNeural = rawDataByLapNeural.rawDataByLapNeural;
                        end
                    elseif strcmp(params.decoding.decID{id}, 'trial')
                        neuralDataFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                            [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'rawDataByTrialNeural.mat');
                        if isfile(neuralDataFname)
                            rawDataByTrialNeural = load(neuralDataFname);
                            rawDataByTrialNeural = rawDataByTrialNeural.rawDataByTrialNeural;
                        end
                    else
                        print('Incorrect name of neural file for decoding.')
                    end
                    %session data for putative cell types
                    sessionDataFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                        [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'rawDataBySessionNeural.mat');
                    if isfile(sessionDataFname)
                        sessionData = load(sessionDataFname);
                        sessionData = sessionData.rawDataBySessionNeural;
                    end
                    %file info for zone info
                    fileInfoFname = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData\' params.iden num2str(animal)], ...
                        [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'fileInfo.mat');
                    if isfile(fileInfoFname)
                        fileInfo = load(fileInfoFname);
                        fileInfo = fileInfo.virmen_fileInfo;
                    end
                    [params.Azones, params.Rzones, params.NRzones, params.AltNRzones] = getZoneInfo_linearJLK(fileInfo, [params.iden num2str(animal)]);

                    if ~isempty(sessionData)
                        %%%%% name the current environment %%%%%
                        if sessionType == 1 %session in original environment
                            currEnv = 'og';
                        elseif sessionType == 2 %session in update environment
                            currEnv = 'up';
                        elseif sessionType == 3 %session in novel environment
                            currEnv = 'nov';
                        elseif sessionType == 4 %session in novel2 environment
                            currEnv = 'nov2';
                        else
                            disp(['Session type not specified for day ' num2str(sessionInfo(f,2))])
                        end

                        %%%%% determine cells to include %%%%%
                        %initialize with no cells included
                        cells2include = zeros(1,length(sessionData.apData));
                        %cell types: all/pyramidal/interneuron
                        if strcmp(params.decoding.cellType, 'all')
                            cellType = ones(1,length(sessionData.apData));%all cells
                        elseif strcmp(params.decoding.cellType, 'pyr')
                            cellType = zeros(1,length(sessionData.apData));
                            for clu = 1:length(sessionData.apData)
                                cellType(clu) = contains(sessionData.apData(clu).putativeCellType, 'Pyr');
                            end
                        elseif strcmp(params.decoding.cellType, 'narrow')
                            cellType = zeros(1,length(sessionData.apData));
                            for clu = 1:length(sessionData.apData)
                                cellType(clu) = contains(sessionData.apData(clu).putativeCellType, 'Narrow');
                            end
                        end%cell type
                        %cell location: CA1/CA3/HIP/cortical
                        if strcmp(params.decoding.cellLocation, 'all')
                            cellLocation = ones(1,length(sessionData.apData));%all cells
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            %%%%%RETURN TO ADD CELL LOCATION ELSEIFs%%%%
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        end%cell location
                        %index by filters
                        cells2include(cellType & cellLocation) = 1;

                        %%%%% collect data %%%%%
                        if strcmp(params.decoding.decID{id}, 'lap')
                            data.(currEnv).apData(ss).ID = [rawDataByLapNeural(1).apData(find(cells2include)).ID];
                            data.(currEnv).apData(ss).maxChan = [rawDataByLapNeural(1).apData(find(cells2include)).maxChan];
                            data.(currEnv).apData(ss).putativeCellType = [rawDataByLapNeural(1).apData(find(cells2include)).putativeCellType];
                            data.(currEnv).degBinEdges = rawDataByLapNeural(1).degBinEdges;
                            for lp = 1:length(rawDataByLapNeural)
                                if length(rawDataByLapNeural(lp).degBinEdges) == length(rawDataByLapNeural(1).degBinEdges)%ensure trial has correct number of bins
                                    lpCtr(1) = lpCtr(1)+1;
                                    for rz = 1:length(params.Rzones)
                                        data.(currEnv).numZoneRewards(ss,lpCtr(1),rz) = sum(rawDataByLapNeural(lp).rewarded(rawDataByLapNeural(lp).currentDeg >= params.Rzones(rz) &...
                                            rawDataByLapNeural(lp).currentDeg <= params.Rzones(rz)+params.cueSize));
                                        isCorrectTmp(rz) = double(data.(currEnv).numZoneRewards(ss,lpCtr(1),rz)>1);%use if params.decoding.correctCriterion = 'reward'
                                    end%rz
                                    %determine correct/incorrect lap
                                    if strcmp(params.decoding.correctCriterion, 'reward')
                                        data.(currEnv).isCorrect(ss,lpCtr(1)) = double(sum(isCorrectTmp)==length(params.Rzones));
                                    elseif strcmp(params.decoding.correctCriterion, 'anticipation')
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        %%%%%ADD ANTICIPATION HERE USING STATSBYREWARDTRIAL%%%%%
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    end%correct/incorrect
                                    data.(currEnv).degBinOccupSmooth(ss,lpCtr(1),:) = rawDataByLapNeural(lp).degBinOccupSmooth;
                                    data.(currEnv).degBinSpeedSmooth(ss,lpCtr(1),:) = rawDataByLapNeural(lp).degBinSpeedSmooth;
                                    data.(currEnv).degBinSpikeCountSmooth(ss,lpCtr(1),1:size([rawDataByLapNeural(lp).degBinSpikeCountSmooth(:,find(cells2include))],1),...
                                        1:size([rawDataByLapNeural(lp).degBinSpikeCountSmooth(:,find(cells2include))],2)) =...
                                        [rawDataByLapNeural(lp).degBinSpikeCountSmooth(:,find(cells2include))];
                                    data.(currEnv).degBinRateMap(ss,lpCtr(1),1:size([rawDataByLapNeural(lp).degBinRateMap(:,find(cells2include))],1),...
                                        1:size([rawDataByLapNeural(lp).degBinRateMap(:,find(cells2include))],2)) =...
                                        [rawDataByLapNeural(lp).degBinRateMap(:,find(cells2include))];
                                end%length of degBinEdges
                            end%lp
                        elseif strcmp(params.decoding.decID{id}, 'trial')
                            data.(currEnv).apData(ss).ID = [rawDataByTrialNeural{1,end}(1).apData(find(cells2include)).ID];
                            data.(currEnv).apData(ss).maxChan = [rawDataByTrialNeural{1,end}(1).apData(find(cells2include)).maxChan];
                            data.(currEnv).apData(ss).putativeCellType = [rawDataByTrialNeural{1,end}(1).apData(find(cells2include)).putativeCellType];
                            data.(currEnv).degBinEdges = rawDataByTrialNeural{1,end}(1).degBinEdges;
                            for znType = 1:size(rawDataByTrialNeural,1)%1=reward, 2=nonreward, 3 = alt nonreward
                                for lp = 1:length(rawDataByTrialNeural{znType,end})
                                    if length(rawDataByTrialNeural{znType,end}(lp).degBinEdges) == length(rawDataByTrialNeural{znType,end}(1).degBinEdges)%ensure trial has correct number of bins
                                        lpCtr(znType) = lpCtr(znType)+1;
                                        if znType == 1
                                            data.(currEnv).numZoneRewards(ss,lpCtr(znType)) = sum(rawDataByTrialNeural{znType,end}(lp).rewarded);
                                            %determine correct/incorrect trial
                                            if strcmp(params.decoding.correctCriterion, 'reward')
                                                data.(currEnv).isCorrect(ss,lpCtr(znType)) = double(data.(currEnv).numZoneRewards(ss,lpCtr(znType))>1);
                                            elseif strcmp(params.decoding.correctCriterion, 'anticipation')
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                %%%%%ADD ANTICIPATION HERE USING STATSBYREWARDTRIAL%%%%%
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            end%correct/incorrect
                                        end%rz
                                        data.(currEnv).degBinOccupSmooth(ss,znType,lpCtr(znType),:) = rawDataByTrialNeural{znType,end}(lp).degBinOccupSmooth;
                                        data.(currEnv).degBinSpeedSmooth(ss,znType,lpCtr(znType),:) = rawDataByTrialNeural{znType,end}(lp).degBinSpeedSmooth;
                                        data.(currEnv).degBinSpikeCountSmooth(ss,znType,lpCtr(znType),1:size([rawDataByTrialNeural{znType,end}(lp).degBinSpikeCountSmooth(:,find(cells2include))],1),...
                                            1:size([rawDataByTrialNeural{znType,end}(lp).degBinSpikeCountSmooth(:,find(cells2include))],2)) =...
                                            [rawDataByTrialNeural{znType,end}(lp).degBinSpikeCountSmooth(:,find(cells2include))];
                                        data.(currEnv).degBinRateMap(ss,znType,lpCtr(znType),1:size([rawDataByTrialNeural{znType,end}(lp).degBinRateMap(:,find(cells2include))],1),...
                                            1:size([rawDataByTrialNeural{znType,end}(lp).degBinRateMap(:,find(cells2include))],2)) =...
                                            [rawDataByTrialNeural{znType,end}(lp).degBinRateMap(:,find(cells2include))];
                                    end%length of degBinEdges
                                end%lp
                            end%%znType
                        end%name of params.decoding.decID

                    end%if ~isempty(neuralData)

                    %save the zone and session info and decoding parameters
                    data.(currEnv).azBins_deg = params.Azones;
                    data.(currEnv).czBins_deg = params.NRzones;
                    data.(currEnv).altczBins_deg = params.AltNRzones;
                    data.(currEnv).sessionInfo{ss} = sessionInfo;
                    data.(currEnv).decodingParams.trialtype = params.decoding.trialType;
                    data.(currEnv).decodingParams.cellType = params.decoding.cellType;
                    data.(currEnv).decodingParams.cellLocation = params.decoding.cellLocation;
                    data.(currEnv).decodingParams.correctCriterion = params.decoding.correctCriterion;

                end%files

                %%%%% determine trials to include %%%%%
                trials2include = zeros(ss,length(data.(currEnv).isCorrect(ss,:)));
                if strcmp(params.decoding.trialType, 'all')
                    trials2include = ones(ss,length(data.(currEnv).isCorrect));
                elseif strcmp(params.decoding.trialType, 'correct')
                    trials2include(data.(currEnv).isCorrect==1) = 1;
                elseif strcmp(params.decoding.trialType, 'incorrect')
                    trials2include(data.(currEnv).isCorrect==0) = 1;
                end
                %filter data
                if strcmp(params.decoding.decID{id}, 'lap')
                    data.(currEnv).degBinOccupSmooth(ss,logical(trials2include(ss,:))==0,:) = 0;
                    data.(currEnv).degBinSpikeCountSmooth(ss,logical(trials2include(ss,:))==0,:,:) = 0;
                    data.(currEnv).degBinRateMap(ss,logical(trials2include(ss,:))==0,:,:) = 0;
                elseif strcmp(params.decoding.decID{id}, 'trial')
                    data.(currEnv).degBinOccupSmooth(ss,:,logical(trials2include(ss,:))==0,:) = 0;
                    data.(currEnv).degBinSpikeCountSmooth(ss,:,logical(trials2include(ss,:))==0,:,:) = 0;
                    data.(currEnv).degBinRateMap(ss,:,logical(trials2include(ss,:))==0,:,:) = 0;
                end%name of params.decoding.decID

            end%~isempty(sessionInfo)

            sprintf('Finished gathering decoding data for %s%d_%d', params.iden, uniqSess(ss,1),uniqSess(ss,2))

        end%session

        %save data
        dir2save = fullfile(dirs.saveoutputstructs, 'Data\Neural\Decoding');
        if ~isfolder(dir2save); mkdir(dir2save); end
        save([dir2save, '\', params.decoding.decID{id}, 'Data.mat'], 'data', '-v7.3');

    else
        load(decDatafname)
    end%~exist(decDatafname) || params.rewrite.decoding.data

end%id


%% perform decoding with data %%

for id = 1%:length(params.decoding.decID)

    decResultsfname = fullfile(dirs.saveoutputstructs, 'Data\Neural\Decoding', [params.decoding.decID{id} 'Results.mat']);

    if ~exist(decResultsfname) || params.rewrite.decodingResults

        %load data
        load([dirs.saveoutputstructs, 'Data\Neural\Decoding\', params.decoding.decID{id} 'Data.mat'])

        %loop through unique sessions
        for ss = 7%1:8%size(uniqSess,1)

            %find number of included trials
            if strcmp(params.decoding.decID{id}, 'lap')
                includedTrials = find(data.(currEnv).degBinOccupSmooth(ss,:,1)~=0);
            elseif strcmp(params.decoding.decID{id}, 'trial')
                includedTrials = find(data.(currEnv).degBinOccupSmooth(ss,1,:,1)~=0);
            end%strcmp

            %loop through environments
            for ee = 1%:length(params.environments)
                currEnv = params.environments{ee};

                decOut = [];
                decOutMat = [];
                if strcmp(params.decoding.analysisType, 'kfold')
                    if length(includedTrials) > params.decoding.minTrials
                        indices = crossvalind('Kfold', length(includedTrials), params.decoding.kfold);
                        for k = 1:params.decoding.kfold
                            %%%%% get training/testing trials %%%%%
                            trainInd = find(indices ~= k);
                            testInd = find(indices == k);
                            trainTrialInd = includedTrials(trainInd);
                            testTrialInd = includedTrials(testInd);

                            if strcmp(params.decoding.decID{id}, 'lap')
                                %get training data
                                trainingData = squeeze(nanmean(data.(currEnv).degBinRateMap(ss,trainTrialInd,:,:),2));
                                trainingData = round(trainingData .* params.decoding.scaleFactor);%make integers
                                tmpTrainingData = [];
                                tmpTrainingData = transpose(squeeze(trainingData(:,:)));%cell x position bin

                                %loop through test trials
                                for tr = 1:length(testTrialInd)
                                    %get testing data
                                    testingData = transpose(squeeze(data.(currEnv).degBinRateMap(ss,testTrialInd(tr),:,:)));%cell x position bin
                                    testingData = round(testingData .* params.decoding.scaleFactor);%make integers

                                    %%%%% DO DECODING %%%%%
                                    %adapted from getspatialprob.m

                                    spatialProb = [];
                                    isLowOcc = isnan(tmpTrainingData(1,:));
                                    for b = 1:size(testingData,2)
                                        binSpikeCount = repmat(testingData(:,b), 1, size(tmpTrainingData,2));
                                        prob = prod(((tmpTrainingData.^binSpikeCount)./gamma(binSpikeCount + 1)).*exp(-tmpTrainingData),1)'; % XZ modified 2/23/2025

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
                                    end
                        
                                    %save to output structs
                                    decOut(k,tr,:).spatialProb = spatialProb;%bins x time
                                    decOut(k,tr,:).totalSpikes = sum(testingData,1);
                                    decOutMat(k,tr,:,:) = spatialProb;
                                end%tr

                            elseif strcmp(params.decoding.decID{id}, 'trial')
                                trainingData = squeeze(nanmean(data.(currEnv).degBinRateMap(ss,:,trainTrialInd,:,:),3));
                                trainingData = round(trainingData .* params.decoding.scaleFactor);%make integers
                                %loop through trial types
                                for znType = 1:size(trainingData,1)
                                    %get training data for this zone type
                                    tmpTrainingData = [];
                                    tmpTrainingData = transpose(squeeze(trainingData(znType,:,:)));

                                    %loop through test trials
                                    for tr = 1:length(testTrialInd)
                                        %get testing data
                                        testingData = transpose(squeeze(data.(currEnv).degBinRateMap(ss,znType,testTrialInd(tr),:,:)));
                                        testingData = round(testingData .* params.decoding.scaleFactor);%make integers

                                        %%%%% DO DECODING %%%%%
                                        %adapted from getspatialprob.m

                                        isLowOcc = isnan(tmpTrainingData(1,:));
                                        for b = 1:size(testingData,2)
                                            binSpikeCount = repmat(testingData(:,b), 1, size(tmpTrainingData,2));
                                            prob = prod(((tmpTrainingData.^binSpikeCount)./gamma(binSpikeCount + 1)).*exp(-tmpTrainingData),1)'; % XZ modified 2/23/2025

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
                                        end

                                        %save to output structs
                                        decOut(k,znType,tr,:).spatialProb = spatialProb; %bins x time
                                        decOut(k,znType,tr,:).totalSpikes = sum(testingData,1);
                                        decOutMat(k,znType,tr,:,:) = spatialProb;

                                    end%tr
                                end%znType
                            end%trials or laps
                        end%k

                    else
                        sprintf('Not enough included %s %ss (<%d) for JK%d_%d, only %d.', currEnv, params.decoding.decID{id}, ...
                            params.decoding.minTrials, data.(currEnv).sessionInfo{ss}(1,1), data.(currEnv).sessionInfo{ss}(1,2), length(includedTrials))
                    end%included trials > minTrials
                elseif strcmp(params.decoding.analysisType, 'leave1out')
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%ADD LEAVE1OUT METHOD%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end%params.decoding.analysisType

            end%ee
        end%ss
    end% ~exist(decResultsfname) || params.rewrite.decodingResults
end%id


%for laps
imagesc(squeeze(nanmean(nanmean(decOutMat,2),1)))
%for trials
imagesc(squeeze(nanmean(nanmean(decOutMat(:,2,:,:,:),3),1)))









%save data
dir2save = fullfile(dirs.saveoutputstructs, 'Data\Behavior\ROC', [params.iden num2str(animal)], ...
    num2str(sessionInfo(1,2)));
if ~isfolder(dir2save); mkdir(dir2save); end
save([dir2save, '\', params.rocID{id}, '.mat'], 'data', 'data');


else
    load(decfname);
end