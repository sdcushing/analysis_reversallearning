function [rateMapAll] = getFiringRateResultsByDay_linearJLK(dirs,uniqSess,params)
%code follows similar structure as getBehaviorROC_JLK.m,
% getSpeedandLickrateResultsByDay_linearJLK, and getDecodingResults_linearJLK

%% Create struct for individual days %%
for id = 1:2%:length(params.decoding.decID)

    %loop through unique sessions
    for ss = 1:51%size(uniqSess,1)

        dayDatafname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' params.decoding.decID{id} 'Data.mat']);

        %loop through environments
        for ee = 1:2%:length(params.environments)
            currEnv = params.environments{ee};

            if isfile(dayDatafname)
                %load this day's neural data
                data = load(dayDatafname);
                data = data.data;

                %only loop if session exists (included because og sessions
                % have fewer entries than update sessions)
                if ~isempty(data.(currEnv))

                    for cellT = 1:2%1=pyr, 2 = narrow spiking
                        %initialize structs for this session
                        rateMap = [];
                        rateMapfname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\' params.iden...
                            num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' currEnv params.decoding.decID{id} params.decoding.cellTypeOptions{cellT} 'RateMap.mat']);

                        if strcmp(params.decoding.decID{id}, 'lap')
                            %loop through cells
                            cluCtr = 0;
                            %determine degree bins
                            tmpDegEdges = data.(currEnv).degBinEdges;
                            tmpDegEdges(1) = 0.0001;
                            for clu = 1:length(data.(currEnv).apData(end).ID)
                                %ONLY DOING NARROW SPIKING INTERNEURONS FOR NOW
                                if contains(data.(currEnv).apData(end).putativeCellType{clu}, params.decoding.cellTypeOptions{cellT}, 'IgnoreCase', true)
                                    cluCtr = cluCtr + 1;
                                    %create ratemap: collect spikes per spatial bin
                                    % across all trials, sum occupancy per spatial bin
                                    % across trials, smooth each, then divide smoothed spike
                                    % counts by smoothed occupancy
                                    tmpSpikeCounts = []; tmpOccup = []; tmpRateMap = [];
                                    tmpSpikeCounts = gaussSmooth(histcounts(squeeze(data.(currEnv).spikePosAll(end,:,:,clu)), tmpDegEdges), 5);
                                    tmpOccup = gaussSmooth(squeeze(nansum(data.(currEnv).degBinOccup(end, :, :), 2))', 5);
                                    tmpRateMap = tmpSpikeCounts ./ tmpOccup;

                                    %regress out smoothed speed and lick rate
                                    ogRateMapSize = size(tmpRateMap);
                                    speedFlat = reshape(squeeze(nanmean(data.(currEnv).degBinSpeedSmooth(end,:,:))),[],1);
                                    lickRateFlat = reshape(squeeze(nanmean(data.(currEnv).degBinLickRateSmooth(end,:,:))),[],1);
                                    tmpRateMapFlat = reshape(tmpRateMap,[],1);
                                    X = [speedFlat, lickRateFlat, speedFlat.*lickRateFlat, ones(size(speedFlat))];
                                    [~, ~, r, ~, stats] = regress(tmpRateMapFlat, X);
                                    rateMap(cluCtr,:) = r';%cell x position bin
                                end%neuron type
                            end%clu

                            %save data
                            save(rateMapfname, 'rateMap', '-v7.3');

                        elseif strcmp(params.decoding.decID{id}, 'trial')
                            %loop through zones
                            if ~isnan(data.(currEnv).nevrzBins_deg(1)); zns = 3; else zns = 2; end
                            for znType = 1:zns
                                currRZDeg = nan(1,zns);
                                if znType == 1
                                    currRZDeg = data.(currEnv).azBins_deg+10;
                                elseif znType == 2
                                    currRZDeg = data.(currEnv).czBins_deg;
                                elseif znType == 3
                                    currRZDeg = data.(currEnv).nevrzBins_deg;
                                end

                                for znNum = 1:length(currRZDeg)
                                    %degrees for this zone+zone number
                                    if floor(currRZDeg(znNum)+params.cueSize+params.gapAfter-params.binsize_deg) > 360%last bin > 360 deg, use wrapTo360
                                        %determine degree bins
                                        degBinEdges = currRZDeg(znNum)-params.gapBefore:params.binsize_deg:360;
                                        degBinEdges = [degBinEdges 0:params.binsize_deg:wrapTo360(currRZDeg(znNum)+params.cueSize+params.gapAfter-params.binsize_deg)];
                                        degBinEdges(degBinEdges==0) = 0.0001;
                                        [degBinEdges, degBinEdgesIntInd] = sort(degBinEdges);
                                        %remove 0.0001 and 360
                                        degBinEdgesIntInd(degBinEdgesIntInd==1) = [];
                                        degBinEdgesIntInd(degBinEdgesIntInd==max(degBinEdgesIntInd)) = [];
                                        degBinEdgesIntInd = degBinEdgesIntInd - 1;
                                    else%last bin <= 360 deg
                                        degBinEdges = floor(currRZDeg(znNum)-params.gapBefore:params.binsize_deg:currRZDeg(znNum)+params.cueSize+params.gapAfter-params.binsize_deg);
                                    end%floor
                                    if isempty(degBinEdges)
                                        error('degBinEdges is empty, stopping execution.');
                                    end%isempty
                                    %trials for this zone+zone number
                                    znTrials = znNum:length(currRZDeg):size(data.(currEnv).spikePosAll,3);

                                    %loop through cells
                                    cluCtr = 0;
                                    for clu = 1:length(data.(currEnv).apData(end).ID)
                                        %ONLY DOING NARROW SPIKING INTERNEURONS FOR NOW
                                        if contains(data.(currEnv).apData(end).putativeCellType{clu}, params.decoding.cellTypeOptions{cellT}, 'IgnoreCase', true)
                                            cluCtr = cluCtr + 1;
                                            %create ratemap: collect spikes per spatial bin
                                            % across all trials, sum occupancy per spatial bin
                                            % across trials, smooth each, then divide smoothed spike
                                            % counts by smoothed occupancy
                                            tmpSpikeCounts = []; tmpOccup = []; tmpRateMap = [];
                                            if floor(currRZDeg(znNum)+params.cueSize+params.gapAfter-params.binsize_deg) > 360%last bin > 360 deg, use wrapTo360
                                                tmpSpikeCounts = histcounts(squeeze(data.(currEnv).spikePosAll(end,znType,znTrials,:,clu)), degBinEdges(degBinEdges<300));
                                                tmpSpikeCounts = [tmpSpikeCounts histcounts(squeeze(data.(currEnv).spikePosAll(end,znType,znTrials,:,clu)), degBinEdges(degBinEdges>300))];
                                                tmpSpikeCounts = tmpSpikeCounts(degBinEdgesIntInd);
                                                tmpSpikeCounts = gaussSmooth(tmpSpikeCounts, 2);
                                            else%last bin <= 360 deg
                                                tmpSpikeCounts = gaussSmooth(histcounts(squeeze(data.(currEnv).spikePosAll(end,znType,znTrials,:,clu)), degBinEdges), 2);
                                            end
                                            tmpOccup = gaussSmooth(squeeze(nansum(data.(currEnv).degBinOccup(end,znType,znTrials,1:length(tmpSpikeCounts)), 3))', 2);
                                            tmpRateMap = tmpSpikeCounts ./ tmpOccup;

                                            %regress out smoothed speed and lick rate
                                            ogRateMapSize = size(tmpRateMap);
                                            speedFlat = reshape(squeeze(nanmean(data.(currEnv).degBinSpeedSmooth(end,znType,znTrials,1:length(tmpSpikeCounts)))),[],1);
                                            lickRateFlat = reshape(squeeze(nanmean(data.(currEnv).degBinLickRateSmooth(end,znType,znTrials,1:length(tmpSpikeCounts)))),[],1);
                                            tmpRateMapFlat = reshape(tmpRateMap,[],1);
                                            X = [speedFlat, lickRateFlat, speedFlat.*lickRateFlat, ones(size(speedFlat))];
                                            [~, ~, r, ~, stats] = regress(tmpRateMapFlat, X);
                                            rateMap(znType,cluCtr,:) = r';%cell x position bin
                                        end%neuron type
                                    end%clu
                                end%znNum
                            end%znType

                            %save data
                            save(rateMapfname, 'rateMap', '-v7.3');

                        end%if lap/trial

                        sprintf('Finished gathering %s firing rates per %s per day for %s%d_%d', params.decoding.cellTypeOptions{cellT}, params.decoding.decID{id}, params.iden, uniqSess(ss,1), uniqSess(ss,2))

                    end%cellT
                end%if ~isempty(data.(currEnv))
            end%if isfile(dayDatafname)
        end%ee
    end%ss
end%id

%% Create struct for all days %%
%done separately from individual days due to computational restraints
for id = 1:2%:length(params.decoding.decID)
    for cellT = 1:2%1=pyr, 2 = narrow spiking
        %initialize output structs
        if strcmp(params.decoding.decID{id}, 'lap')
            %output = sessions x environments x cell x position bin
            rateMapAll = nan(size(uniqSess,1),2,500,180);
        elseif strcmp(params.decoding.decID{id}, 'trial')
            %output = sessions x environments x zone type x cell x position bin
            rateMapAll = nan(size(uniqSess,1),2,3,500,50);
        end%if lap/trial

        %loop through unique sessions
        for ss = 1:51%size(uniqSess,1)

            %loop through environments
            for ee = 1:2%:length(params.environments)
                currEnv = params.environments{ee};

                rateMapfname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' currEnv params.decoding.decID{id}  params.decoding.cellTypeOptions{cellT} 'RateMap.mat']);

                if isfile(rateMapfname)

                    %load this day's data
                    rateMap = [];
                    rateMap = load(rateMapfname);
                    rateMap = rateMap.rateMap;

                    %only loop if session exists (included because og sessions
                    % have fewer entries than update sessions)
                    if ~isempty(rateMap)

                        if strcmp(params.decoding.decID{id}, 'lap')
                            %add data to struct for all sessions
                            % output = sessions x environments x cell x position bin
                            rateMapAll(ss,ee,1:size(rateMap,1),:) = rateMap;
                        elseif strcmp(params.decoding.decID{id}, 'trial')
                            %add data to struct for all sessions
                            % output = sessions x environments x zone type x cell x position bin
                            rateMapAll(ss,ee,1:size(rateMap,1),1:size(rateMap,2),1:size(rateMap,3)) = rateMap;
                        end%if lap/trial

                    end%if ~isempty(data.(currEnv))
                end%if isfile(rateMapfname)
            end%ee

        end%ss

        %save data
        rateMapAllfname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\' params.decoding.decID{id}  params.decoding.cellTypeOptions{cellT} 'RateMapAll.mat']);
        save(rateMapAllfname, 'rateMapAll', '-v7.3');

        sprintf('Finished gathering %s firing rates per %s for all days',  params.decoding.cellTypeOptions{cellT}, params.decoding.decID{id})

    end%cellT
end%id

end%function