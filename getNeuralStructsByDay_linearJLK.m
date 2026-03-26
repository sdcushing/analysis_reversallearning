function [data] = getNeuralStructsByDay_linearJLK(allindex,dirs,uniqSess,params)
%created by JLK 260203 to collect neural data across sessions on a given
% day; data is filterd by the desired cell type/location, trial type
% (correct/incorrect/all), etc.

%% create data across sessions per day %%

for id = 2%:length(params.decoding.decID)

    %loop through unique sessions
    for ss = 12%1:51%length(uniqSess)

        animal = uniqSess(ss, 1);
        sessDate = uniqSess(ss, 2);
        if id == 1 || id == 2
            tempidx = find(ismember(allindex(:,1), animal) & ismember(allindex(:,2), sessDate) & ismember(allindex(:,6), 2) & ismember(allindex(:,5), 1:4));%include this animal, this date, active and recording only
        elseif id == 3
            tempidx = find(ismember(allindex(:,1), animal) & ismember(allindex(:,2), sessDate) & ismember(allindex(:,6), 2:3) & ismember(allindex(:,5), 1:4));%include this animal, this date, active OR rest and recording only
        end
        sessionInfo = allindex(tempidx, :);

        if ~isempty(sessionInfo)

            %output info
            dayDataDir = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\JK' num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2))]);
            dayDatafname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\JK' num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' params.decoding.decID{id} 'Data.mat']);

            if ~exist(dayDatafname) || params.rewrite.decodingData

                %%%%% make save path if not already created %%%%%
                if ~isfolder(dayDataDir)
                    mkdir(dayDataDir)
                end

                %%% initialize data structure with empty fields %%%%%
                data = [];
                for ee = 1:length(params.environments)
                    data.(params.environments{ee}) = [];
                end

                %loop through files to COLLECT DATA ACROSS SESSIONS FOR THIS DAY
                lpCtr = zeros(3,1);
                lpRipCtr = zeros(4,150,20);
                files = sessionInfo(:,3);
                for f = 1:length(files)
                    currFile = files(f);
                    sessionType = sessionInfo(f,4);
                    %for update days, reset the lap counters for the RZ and NRZ
                    if sum(diff(sessionInfo(:,4)))~=0 && f == find(diff(sessionInfo(:,4)))+1
                        lpCtr = zeros(3,1); lpRipCtr = zeros(4,100,20);
                    end%update day counter reset

                    %%%%% load data %%%%%
                    %lap or trial data for decoding
                    if strcmp(params.decoding.decID{id}, 'lap')
                        neuralDataFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                            [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'rawDataByLapNeural.mat');
                        if isfile(neuralDataFname)
                            rawDataByLapNeural = load(neuralDataFname);
                            rawDataByLapNeural = rawDataByLapNeural.rawDataByLapNeural;
                            tmpCheckData = rawDataByLapNeural;
                        else
                            sprintf('No lap neural file for decoding.')
                        end
                    elseif strcmp(params.decoding.decID{id}, 'trial')
                        neuralDataFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                            [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'rawDataByTrialNeural.mat');
                        if isfile(neuralDataFname)
                            rawDataByTrialNeural = load(neuralDataFname);
                            rawDataByTrialNeural = rawDataByTrialNeural.rawDataByTrialNeural;
                            tmpCheckData = rawDataByTrialNeural;
                        else
                            sprintf('No trial neural file for decoding.')
                        end
                    elseif strcmp(params.decoding.decID{id}, 'ripple') && f < length(files)
                        neuralDataFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                            [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'rawDataByTrialNeural.mat');
                        if isfile(neuralDataFname)
                            rawDataByTrialNeural = load(neuralDataFname);
                            rawDataByTrialNeural = rawDataByTrialNeural.rawDataByTrialNeural;
                            tmpCheckData = rawDataByTrialNeural;
                        else
                            sprintf('No trial neural file for ripple decoding.')
                        end
                    elseif strcmp(params.decoding.decID{id}, 'ripple') && f == length(files)
                        %only need session data below
                    else
                        sprintf('Incorrect name of neural file for decoding.')
                    end%load lap/trial data
                    %session data for putative cell types
                    sessionDataFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                        [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'rawDataBySessionNeural.mat');
                    if isfile(sessionDataFname)
                        sessionData = load(sessionDataFname);
                        sessionData = sessionData.rawDataBySessionNeural;
                        if strcmp(params.decoding.decID{id}, 'ripple') && f == length(files)%rest session only includes session neural struct
                            tmpCheckData = rawDataByTrialNeural;
                        end
                    end%load session data
                    %pyramidal layer info
                    pyrLayersFname = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' params.iden num2str(animal)], ...
                        [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'sessionPyrLayerInfo.mat');
                    if isfile(pyrLayersFname)
                        sessionPyrLayerInfo = load(pyrLayersFname);
                        sessionPyrLayerInfo = sessionPyrLayerInfo.sessionPyrLayerInfo;
                    end%load pyramidal layer
                    %file info for zone info
                    fileInfoFname = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData\' params.iden num2str(animal)], ...
                        [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'fileInfo.mat');
                    if isfile(fileInfoFname)
                        fileInfo = load(fileInfoFname);
                        fileInfo = fileInfo.virmen_fileInfo;
                    end%load zone info
                    [params.Azones, params.Rzones, params.NRzones, params.NevRzones] = getZoneInfo_linearJLK(fileInfo, [params.iden num2str(animal)]);

                    if ~isempty(tmpCheckData)
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
                        cellType = zeros(1,length(sessionData.apData));
                        cellLocation = zeros(1,length(sessionData.apData));
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
                        elseif strcmp(params.decoding.cellLocation, 'HIP')
                            for clu = 1:length(sessionData.apData)
                                cellLocation(clu) = sessionData.apData(clu).maxChan<150;%chan 150 usually roughly top of HIP
                            end
                        elseif strcmp(params.decoding.cellLocation, 'CA1')
                            for clu = 1:length(sessionData.apData)
                                cellLocation(clu) = sessionPyrLayerInfo.pyrLayerCA1(1)<=sessionData.apData(clu).maxChan &...
                                    sessionData.apData(clu).maxChan<=sessionPyrLayerInfo.pyrLayerCA1(end);
                            end
                        elseif strcmp(params.decoding.cellLocation, 'CA3')
                            for clu = 1:length(sessionData.apData)
                                cellLocation(clu) = sessionPyrLayerInfo.pyrLayerCA3(1)<=sessionData.apData(clu).maxChan &...
                                    sessionData.apData(clu).maxChan<=sessionPyrLayerInfo.pyrLayerCA3(end);
                            end
                        end%cell location
                        %index by filters
                        cells2include(cellType & cellLocation) = 1;

                        %%%%% collect data %%%%%
                        %determine start lap/trial if session started recording late
                        stLap = 1;
                        stTr = 1;
                        if sessionInfo(1,1) == 15
                            if sessionInfo(2,2) == 250529 && sessionInfo(2,3) == 3
                                stLap = 4;
                                stTr = 10;
                            end
                        end%stTr JK23_251024_2
                        if sessionInfo(1,1) == 23
                            if sessionInfo(1,2) == 251024 && sessionInfo(1,3) == 2
                                stTr = 2;
                            end
                        end%stTr JK23_251024_2

                        %loop through laps or trials
                        if strcmp(params.decoding.decID{id}, 'lap')
                            data.(currEnv).apData(1).ID = [rawDataByLapNeural(stLap).apData(find(cells2include)).ID];
                            data.(currEnv).apData(1).maxChan = [rawDataByLapNeural(stLap).apData(find(cells2include)).maxChan];
                            data.(currEnv).apData(1).putativeCellType = {rawDataByLapNeural(stLap).apData(find(cells2include)).putativeCellType};
                            data.(currEnv).degBinEdges = rawDataByLapNeural(stLap).degBinEdges;
                            for lp = 1:length(rawDataByLapNeural)
                                if length(rawDataByLapNeural(lp).degBinEdges) == length(rawDataByLapNeural(stLap).degBinEdges)%ensure trial has correct number of bins
                                    lpCtr(1) = lpCtr(1)+1;
                                    data.(currEnv).tBinEdges(1,lpCtr(1),1:size(rawDataByLapNeural(lp).tBinEdges,2)) = rawDataByLapNeural(lp).tBinEdges;
                                    data.(currEnv).tBinSpeed(1,lpCtr(1),1:size(rawDataByLapNeural(lp).tBinSpeed,2)) = rawDataByLapNeural(lp).tBinSpeed;
                                    data.(currEnv).tBinActualPos(1,lpCtr(1),1:size(rawDataByLapNeural(lp).tBinActualPos,2)) = rawDataByLapNeural(lp).tBinActualPos;
                                    for rz = 1:length(params.Rzones)
                                        data.(currEnv).numZoneRewards(1,lpCtr(1),rz) = sum(rawDataByLapNeural(lp).rewarded(rawDataByLapNeural(lp).currentDeg >= params.Rzones(rz) &...
                                            rawDataByLapNeural(lp).currentDeg <= params.Rzones(rz)+params.cueSize));
                                        isCorrectTmp(rz) = double(data.(currEnv).numZoneRewards(1,lpCtr(1),rz)>1);%use if params.decoding.correctCriterion = 'reward'
                                    end%rz
                                    %determine correct/incorrect lap
                                    if strcmp(params.decoding.correctCriterion, 'reward')
                                        data.(currEnv).isCorrect(1,lpCtr(1)) = double(sum(isCorrectTmp)==length(params.Rzones));
                                    elseif strcmp(params.decoding.correctCriterion, 'anticipation')
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        %%%%%ADD ANTICIPATION HERE USING STATSBYREWARDTRIAL%%%%%
                                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                    end%correct/incorrect
                                    data.(currEnv).spikePosMov(1,lpCtr(1),1:size(rawDataByLapNeural(lp).spikePosMov,1),1:size(rawDataByLapNeural(lp).spikePosMov,2)) = rawDataByLapNeural(lp).spikePosMov;
                                    data.(currEnv).spikeTimeMov(1,lpCtr(1),1:size(rawDataByLapNeural(lp).spikeTimeMov,1),1:size(rawDataByLapNeural(lp).spikeTimeMov,2)) = rawDataByLapNeural(lp).spikeTimeMov;
                                    data.(currEnv).spikePosAll(1,lpCtr(1),1:size(rawDataByLapNeural(lp).spikePosAll,1),1:size(rawDataByLapNeural(lp).spikePosAll,2)) = rawDataByLapNeural(lp).spikePosAll;
                                    data.(currEnv).spikeTimeAll(1,lpCtr(1),1:size(rawDataByLapNeural(lp).spikeTimeAll,1),1:size(rawDataByLapNeural(lp).spikeTimeAll,2)) = rawDataByLapNeural(lp).spikeTimeAll;
                                    data.(currEnv).degBinOccup(1,lpCtr(1),:) = rawDataByLapNeural(lp).degBinOccup;
                                    data.(currEnv).degBinLickRateSmooth(1,lpCtr(1),:) = rawDataByLapNeural(lp).degBinLickRateSmooth;
                                    data.(currEnv).degBinSpeedSmooth(1,lpCtr(1),:) = rawDataByLapNeural(lp).degBinSpeedSmooth;
                                end%length of degBinEdges
                            end%lp
                        elseif strcmp(params.decoding.decID{id}, 'trial')
                            data.(currEnv).apData(1).ID = [rawDataByTrialNeural{1,end}(stTr).apData(find(cells2include)).ID];
                            data.(currEnv).apData(1).maxChan = [rawDataByTrialNeural{1,end}(stTr).apData(find(cells2include)).maxChan];
                            data.(currEnv).apData(1).putativeCellType = {rawDataByTrialNeural{1,end}(stTr).apData(find(cells2include)).putativeCellType};
                            for znType = 1:size(rawDataByTrialNeural,1)%1=reward, 2=nonreward, 3 = alt nonreward
                                for lp = 1:length(rawDataByTrialNeural{znType,end})
                                    if ~isempty(rawDataByTrialNeural{znType,end}(lp).degBinOccup) && length(rawDataByTrialNeural{znType,end}(lp).degBinEdges) == length(rawDataByTrialNeural{1,end}(stTr+1).degBinEdges)%ensure trial has data and correct number of bins
                                        lpCtr(znType) = lpCtr(znType)+1;
                                        data.(currEnv).degBinEdges(1,znType,lpCtr(znType),:) = rawDataByTrialNeural{znType,end}(lp).degBinEdges;
                                        data.(currEnv).tBinEdges(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).tBinEdges,2)) = rawDataByTrialNeural{znType,end}(lp).tBinEdges;
                                        data.(currEnv).tBinActualPos(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).tBinActualPos,2)) = rawDataByTrialNeural{znType,end}(lp).tBinActualPos;
                                        data.(currEnv).tBinSpeed(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).tBinSpeed,2)) = rawDataByTrialNeural{znType,end}(lp).tBinSpeed;
                                        if znType == 1
                                            data.(currEnv).numZoneRewards(1,lpCtr(znType)) = sum(rawDataByTrialNeural{znType,end}(lp).rewarded);
                                            %determine correct/incorrect trial
                                            if strcmp(params.decoding.correctCriterion, 'reward')
                                                data.(currEnv).isCorrect(1,lpCtr(znType)) = double(data.(currEnv).numZoneRewards(1,lpCtr(znType))>1);
                                            elseif strcmp(params.decoding.correctCriterion, 'anticipation')
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                %%%%%ADD ANTICIPATION HERE USING STATSBYREWARDTRIAL%%%%%
                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                            end%correct/incorrect
                                        end%if znType == 1
                                        data.(currEnv).spikePosMov(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).spikePosMov,1),1:size(rawDataByTrialNeural{znType,end}(lp).spikePosMov,2)) = rawDataByTrialNeural{znType,end}(lp).spikePosMov;
                                        data.(currEnv).spikeTimeMov(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).spikeTimeMov,1),1:size(rawDataByTrialNeural{znType,end}(lp).spikeTimeMov,2)) = rawDataByTrialNeural{znType,end}(lp).spikeTimeMov;
                                        data.(currEnv).spikePosAll(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).spikePosAll,1),1:size(rawDataByTrialNeural{znType,end}(lp).spikePosAll,2)) = rawDataByTrialNeural{znType,end}(lp).spikePosAll;
                                        data.(currEnv).spikeTimeAll(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).spikeTimeAll,1),1:size(rawDataByTrialNeural{znType,end}(lp).spikeTimeAll,2)) = rawDataByTrialNeural{znType,end}(lp).spikeTimeAll;
                                        data.(currEnv).degBinOccup(1,znType,lpCtr(znType),:) = rawDataByTrialNeural{znType,end}(lp).degBinOccup;
                                        data.(currEnv).degBinLickRateSmooth(1,znType,lpCtr(znType),:) = rawDataByTrialNeural{znType,end}(lp).degBinLickRateSmooth;
                                        data.(currEnv).degBinSpeedSmooth(1,znType,lpCtr(znType),:) = rawDataByTrialNeural{znType,end}(lp).degBinSpeedSmooth;
                                    end%length of degBinEdges
                                end%lp
                            end%%znType
                        elseif strcmp(params.decoding.decID{id}, 'ripple')
                            if f < length(files)%non-rest sessions
                                data.(currEnv).apData(1).ID = [rawDataByTrialNeural{1,end}(stTr).apData(find(cells2include)).ID];
                                data.(currEnv).apData(1).maxChan = [rawDataByTrialNeural{1,end}(stTr).apData(find(cells2include)).maxChan];
                                data.(currEnv).apData(1).putativeCellType = {rawDataByTrialNeural{1,end}(stTr).apData(find(cells2include)).putativeCellType};
                                for znType = 1:size(rawDataByTrialNeural,1)%1=reward, 2=nonreward, 3 = alt nonreward
                                    for lp = 1:length(rawDataByTrialNeural{znType,end})
                                        if ~isempty(rawDataByTrialNeural{znType,end}(lp).degBinOccup) && length(rawDataByTrialNeural{znType,end}(lp).degBinEdges) == length(rawDataByTrialNeural{1,end}(stTr+1).degBinEdges)%ensure trial has data and correct number of bins
                                            lpCtr(znType) = lpCtr(znType)+1;
                                            data.(currEnv).degBinEdges(1,znType,lpCtr(znType),:) = rawDataByTrialNeural{znType,end}(lp).degBinEdges;
                                            data.(currEnv).tBinEdges(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).tBinEdges,2)) = rawDataByTrialNeural{znType,end}(lp).tBinEdges;
                                            data.(currEnv).tBinActualPos(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).tBinActualPos,2)) = rawDataByTrialNeural{znType,end}(lp).tBinActualPos;
                                            data.(currEnv).tBinSpeed(1,znType,lpCtr(znType),1:size(rawDataByTrialNeural{znType,end}(lp).tBinSpeed,2)) = rawDataByTrialNeural{znType,end}(lp).tBinSpeed;
                                            if znType == 1
                                                data.(currEnv).numZoneRewards(1,lpCtr(znType)) = sum(rawDataByTrialNeural{znType,end}(lp).rewarded);
                                                %determine correct/incorrect trial
                                                if strcmp(params.decoding.correctCriterion, 'reward')
                                                    data.(currEnv).isCorrect(1,lpCtr(znType)) = double(data.(currEnv).numZoneRewards(1,lpCtr(znType))>1);
                                                elseif strcmp(params.decoding.correctCriterion, 'anticipation')
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                    %%%%%ADD ANTICIPATION HERE USING STATSBYREWARDTRIAL%%%%%
                                                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                end%correct/incorrect
                                            end%rz
                                            for ch = 1:length(rawDataByTrialNeural{znType,end}(lp).ripplesGood)
                                                for r = 1:length(rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).startind)
                                                    lpRipCtr(znType,lp,ch) = lpRipCtr(znType,lp,ch)+1;
                                                    data.(currEnv).midInd(1,znType,lpCtr(znType),ch,lpRipCtr(znType,lp,ch)) = ...
                                                        rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).midind(r);
                                                    data.(currEnv).currentDeg(1,znType,lpCtr(znType),ch,lpRipCtr(znType,lp,ch)) = ...
                                                        rawDataByTrialNeural{znType,end}(lp).currentDeg(find(rawDataByTrialNeural{znType,end}(lp).lfpTime>...
                                                        data.(currEnv).midInd(1,znType,lpCtr(znType),ch,lpRipCtr(znType,lp,ch)),1,'first'));
                                                    data.(currEnv).degBinOccup(1,znType,lpCtr(znType),ch,lpRipCtr(znType,lp,ch)) = ...
                                                        rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).endtime(r)-rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).starttime(r);
                                                    data.(currEnv).ripSpikes(1,znType,lpCtr(znType),ch,lpRipCtr(znType,lp,ch),...
                                                        1:size(rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).ripSpikes(r,:,:),2),1:size(rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).ripSpikes(r,:,:),3))...
                                                        = squeeze(rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).ripSpikes(r,:,:));
                                                end%r
                                            end%ch
                                        end%~isempty lfpTime
                                    end%lp
                                end%znType
                            elseif f == length(files)%rest sessions
                                znType = 4;%save rest as 'zone' as zone 4
                                for ch = 1:length(rawDataBySessionNeural.ripplesGood)
                                    for r = 1:length(rawDataBySessionNeural.ripplesGood(ch).startind)
                                        lpRipCtr(znType,1,ch) = lpRipCtr(znType,1,ch)+1;
                                        data.(currEnv).midInd(1,znType,1,ch,lpRipCtr(znType,1,ch)) = ...
                                            rawDataBySessionNeural.ripplesGood(ch).midind(r);
                                        data.(currEnv).currentDeg(1,znType,1,ch,lpRipCtr(znType,1,ch)) = ...
                                            rawDataBySessionNeural.currentDeg(find(rawDataBySessionNeural.lfpTime>...
                                            data.(currEnv).midInd(1,znType,1,ch,lpRipCtr(znType,1,ch)),1,'first'));
                                        data.(currEnv).degBinOccup(1,znType,1,ch,lpRipCtr(znType,1,ch)) = ...
                                            rawDataBySessionNeural.ripplesGood(ch).endtime(r)-rawDataBySessionNeural.ripplesGood(ch).starttime(r);
                                        data.(currEnv).ripSpikes(1,znType,1,ch,lpRipCtr(znType,1,ch),...
                                            1:size(rawDataBySessionNeural.ripplesGood(ch).ripSpikes(r,:,:),2),1:size(rawDataBySessionNeural.ripplesGood(ch).ripSpikes(r,:,:),3))...
                                            = squeeze(rawDataBySessionNeural.ripplesGood(ch).ripSpikes(r,:,:));
                                    end%r
                                end%ch
                            end%f
                        end%name of params.decoding.decID
                    end%if ~isempty(tmpCheckData)

                    %save the zone and session info and decoding parameters
                    data.(currEnv).azBins_deg = params.Azones;
                    data.(currEnv).czBins_deg = params.NRzones;
                    data.(currEnv).nevrzBins_deg = params.NevRzones;
                    data.(currEnv).sessionInfo{1} = sessionInfo;
                    data.(currEnv).decodingParams.trialtype = params.decoding.trialType;
                    data.(currEnv).decodingParams.cellType = params.decoding.cellType;
                    data.(currEnv).decodingParams.cellLocation = params.decoding.cellLocation;
                    data.(currEnv).decodingParams.correctCriterion = params.decoding.correctCriterion;

                end%files

                %%%%% determine trials to include %%%%%
                trials2include = zeros(1,length(data.(currEnv).isCorrect(1,:)));
                if strcmp(params.decoding.trialType, 'all')
                    trials2include = ones(1,length(data.(currEnv).isCorrect));
                elseif strcmp(params.decoding.trialType, 'correct')
                    trials2include(data.(currEnv).isCorrect==1) = 1;
                elseif strcmp(params.decoding.trialType, 'incorrect')
                    trials2include(data.(currEnv).isCorrect==0) = 1;
                end
                %filter data
                if strcmp(params.decoding.decID{id}, 'lap')
                    data.(currEnv).degBinOccup(1,logical(trials2include(1,:))==0,:) = 0;
                elseif strcmp(params.decoding.decID{id}, 'trial')
                    data.(currEnv).degBinOccup(1,:,logical(trials2include(1,:))==0,:) = 0;
                elseif strcmp(params.decoding.decID{id}, 'ripple')
                    if f < length(files)%non-rest sessions
                        data.(currEnv).degBinOccup(1,:,logical(trials2include(1,:))==0,:) = 0;
                    elseif f == length(files)%rest sessions
                        %keep all 'trials' for rest sessions
                    end
                end%name of params.decoding.decID

            end%~exist(decDatafname) || params.rewrite.decoding.data

            sprintf('Finished gathering neural data per %s per day for %s%d_%d', params.decoding.decID{id}, params.iden, uniqSess(ss,1),uniqSess(ss,2))

            %save data
            save(dayDatafname, 'data', '-v7.3');

        end%if ~isempty(sessInfo)
    end%ss

end%id

end%function