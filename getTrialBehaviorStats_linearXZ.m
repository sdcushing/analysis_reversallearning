function [rawDataByTrial, statsByRewardTrial, statsByRewardSession] = getTrialBehaviorStats_linearXZ(rawDataBySession, virmen_fileInfo, params, saveBehaviorPath)
%Note: only completed laps
%Output: rawDataByTrial{n,m}(t), statsByRewardTrial{n,m}, and statsByRewardSession{n,1}
% where n = zone type (1 = reward, 2 = nonreward, 3 = alt nonreward),
% m = zone # (1,2,3,4), and t = trial #.

%% Get stats by trial %%
%adapted from getTrialByTrialStats_linearJLK

% split data up into individual laps (0 to 360 deg)
%lapStartIdx = index of starting a new lap, going from 0 to 360 deg
lapStartIdx = find(diff(rawDataBySession.currentDeg)<-300)+1;
lapStartIdx = [1;lapStartIdx(1:end)];

% Calculate trial by trial stats
clear numRewards totalTime lickBehavior lickCounts lickProbability lickRate lickRateSmooth velocAvg velocCounts velocCountsSmooth rawDataByTrial statsByRewardTrial statsByRewardSession

if length(lapStartIdx) < 2 %no completed laps

    %save to statsByRewardTrial struct
    statsByRewardTrial = [];

else %at least one completed lap

    %%%%% create raw data struct %%%%%

    %parameters, same for each zone type
    endBin = abs(params.gapBefore) + abs(params.gapAfter);
    binEdges = -params.gapBefore:params.binsize_deg:params.gapAfter;
    histEdges = 0:params.binsize_deg:endBin;

    for znType = 1:3 %1=reward, 2=nonreward, 3 = alt nonreward

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% determine zones %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if znType == 1
            znStartDegs = params.Rzones;
        elseif znType == 2
            znStartDegs = params.NRzones;
        elseif znType == 3
            znStartDegs = params.AltNRzones;
        end

        %padded zones are X degrees (as defined by gap) before and after each reward zone
        % #% XZ DEBUG: COMPATIBLE WITH EARLY/LATE RZs
        paddedZones = {};
        for iZone = 1:length(znStartDegs)   
            paddedZones{iZone} = [wrapTo360(znStartDegs(iZone) - params.gapBefore), wrapTo360(znStartDegs(iZone) + virmen_fileInfo.cueSize + params.gapAfter)];
            % deal with edge case
            if paddedZones{iZone}(1) > paddedZones{iZone}(2)
                seg_1 = [paddedZones{iZone}(1), 360];
                seg_2 = [0, paddedZones{iZone}(2)];
                paddedZones{iZone} = [seg_1; seg_2];
            end
        end
        % #
        % 
        % paddedZones = [];
        % paddedZones(:,1) = wrapTo360(znStartDegs - params.gapBefore);
        % paddedZones(:,2) = wrapTo360(znStartDegs + virmen_fileInfo.cueSize + params.gapAfter);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% split data up into individual trials %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        allC = [];
        for p = 1:length(znStartDegs) %in the order of zones
            % temp = isExcluded(rawDataBySession.currentDeg, paddedZones(p,:));
            temp = isExcluded(rawDataBySession.currentDeg, paddedZones{p});
            ftemp = find(temp);
            if ~isempty(ftemp)
                fidx = find(diff(ftemp)>1);
                fidx = [1; fidx];
                PZ = [ftemp(1);ftemp(fidx)];
                PZ = [PZ; ftemp(fidx+1); ftemp(end)];
                PZ = sort(PZ,'ascend');
                if mod(size(PZ,1),2) == 0 %ensure it's even number of elements
                    C = [PZ(1:2:end),PZ(2:2:end)];
                    %exclude trials that are shorter (in degrees) than expected
                    idx2excl = find(wrapTo360(rawDataBySession.currentDeg(C(:,2)) - rawDataBySession.currentDeg(C(:,1))) < params.gapBefore + params.gapAfter);
                    C(idx2excl,:) = [];
                    allC = [allC; C];
                else
                    C = [PZ(1:2:end-1),PZ(2:2:end-1)];
                    idx2excl = find(wrapTo360(rawDataBySession.currentDeg(C(:,2)) - rawDataBySession.currentDeg(C(:,1))) < params.gapBefore + params.gapAfter);
                    if ~isempty(idx2excl)
                        C(idx2excl,:) = [];
                    end
                    allC = [allC; C];
                end
                for ii = 1:size(C,1)
                    rawDataByTrial{znType,p}(ii) = structfun(@(x) x(C(ii,1):C(ii,2)), rawDataBySession, 'UniformOutput', false);
                end
                clear PZ
            end
            allC = sort(allC,1); %last struct uses all zones in order of appaearance

            for ii = 1:size(allC,1)
                rawDataByTrial{znType,length(paddedZones)+1}(ii) = structfun(@(x) x(allC(ii,1):allC(ii,2)), rawDataBySession, 'UniformOutput', false);
            end

        end%p
    end%znType

    %%%%% stats for trials %%%%%

    %initialize output struct
    statsByRewardTrial = cell(size(rawDataByTrial));

    for znType = 1:size(rawDataByTrial,1)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% determine zones %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if znType == 1%reward
            znStartDegs = params.Rzones;
        elseif znType == 2%control
            znStartDegs = params.NRzones;
        elseif znType == 3%alt control
            znStartDegs = params.AltNRzones;
        end

        for znNum = 1:length(znStartDegs)
            if ~isempty(rawDataByTrial{znType,znNum})

                %initialize variables
                lickCounts = zeros(length(rawDataByTrial{znType,znNum}),length(histEdges)-1);
                lickProbability = zeros(length(rawDataByTrial{znType,znNum}),length(histEdges)-1);
                lickRate = zeros(length(rawDataByTrial{znType,znNum}),length(histEdges)-1);
                velocCounts = zeros(length(rawDataByTrial{znType,znNum}),length(histEdges)-1);
                velocCountsSmooth = zeros(length(rawDataByTrial{znType,znNum}),length(histEdges)-1);

                for tr = 1:length(rawDataByTrial{znType,znNum})

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% some initial variables for the trial %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    lickTimes = [0; diff(rawDataByTrial{znType,znNum}(tr).licks)];

                    timeout = find(rawDataByTrial{znType,znNum}(tr).currentZone == params.timeoutZone);
                    %find times during timeout periods to exclude
                    if ~isempty(timeout)
                        lickTimes(timeout) = 0; %converts any licks that occurred during timeout to zero
                    end
                    totalLicks = sum(lickTimes); %total lickcount per trial

                    trDeg = rawDataByTrial{znType,znNum}(tr).currentDeg;
                    trDeg = [0; wrapTo360(diff(trDeg))];
                    trDeg(trDeg>300) = 0;
                    trDegCum = cumsum(trDeg); %final vel = degrees that go from zero to length of trial window

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% reward, time values %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    numRewards(tr) = max(rawDataByTrial{znType,znNum}(tr).rewards)-min(rawDataByTrial{znType,znNum}(tr).rewards);
                    totalTime(tr) = rawDataByTrial{znType,znNum}(tr).vrTime(end)-rawDataByTrial{znType,znNum}(tr).vrTime(1);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% licking behavior %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %adapted from getLickInfoForRewardTrial_linearJLK

                    %number of licks in each zone
                    lickBehavior(tr).licksInZone = sum(lickTimes(find(rawDataByTrial{znType,znNum}(tr).currentDeg >= znStartDegs(znNum) & rawDataByTrial{znType,znNum}(tr).currentDeg < wrapTo360(znStartDegs(znNum) + virmen_fileInfo.cueSize))));
                    lickBehavior(tr).timeInZone = sum(diff(rawDataByTrial{znType,znNum}(tr).vrTime(find(rawDataByTrial{znType,znNum}(tr).currentDeg >= znStartDegs(znNum) & rawDataByTrial{znType,znNum}(tr).currentDeg < wrapTo360(znStartDegs(znNum) + virmen_fileInfo.cueSize )))));

                    %lick rate per zone
                    lickBehavior(tr).lickRateInZone = lickBehavior(tr).licksInZone/lickBehavior(tr).timeInZone;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% licking counts %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%
                    %adapted from lickCountsForRewardTrial_linearJLK

                    %calculate lick bins across the trial
                    for i = 1:length(histEdges)-1
                        binIdx = find(trDegCum >= histEdges(i) & trDegCum < histEdges(i+1));
                        if ~isempty(binIdx)
                            lickCounts(tr,i) = sum(lickTimes(binIdx));
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% licking probability %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %adapted from lickProbabilityForRewardTrial_linearJLK

                    %calculate lick bins across the trial
                    for i = 1:length(histEdges)-1
                        binIdx = find(trDegCum >= histEdges(i) & trDegCum < histEdges(i+1));
                        if ~isempty(binIdx)
                            lickProbability(tr,i) = sum(lickTimes(binIdx))/totalLicks;
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%
                    %%%%% lick rate %%%%%
                    %%%%%%%%%%%%%%%%%%%%%
                    %adapted from lickRateForTrial_linearJLK

                    %calculate lick bins across the trial
                    for i = 1:length(histEdges)-1
                        binIdx = find(trDegCum >= histEdges(i) & trDegCum < histEdges(i+1));
                        if ~isempty(binIdx)
                            binTime = rawDataByTrial{znType,znNum}(tr).vrTime(binIdx);
                            lickRate(tr,i) = sum(lickTimes(binIdx)) / (binTime(end)-binTime(1));
                        else
                            lickRate(tr,i) = 0;
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% velocity average %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    velocAvg(tr,:) = sum(wrapTo360(diff(rawDataByTrial{znType,znNum}(tr).currentDeg)))/totalTime(tr);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% velocity counts %%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %calculate instantaneous velocity from smoothed coordinates across trial
                    velocVect = zeros(1,length(rawDataByTrial{znType,znNum}(tr).currentDeg));
                    timeDiff = [0; diff(rawDataByTrial{znType,znNum}(tr).vrTime)];

                    for i = 1:length(rawDataByTrial{znType,znNum}(tr).currentDeg)
                        if rawDataByTrial{znType,znNum}(tr).currentZone(i) == params.timeoutZone
                            % exclude time-out periods from velocity calculation
                            velocVect(i) = 0; %timeout velocity printed as zero (no actual movement through VR during timeout period)
                        else
                            velocVect(i) = trDeg(i) ./ timeDiff(i);
                        end
                    end

                    %get velocity bins
                    for i = 1:length(histEdges)-1
                        binIdx = find(trDegCum >= histEdges(i) & trDegCum < histEdges(i+1));
                        if ~isempty(binIdx)
                            velocCounts(tr,i) = nanmean(velocVect(binIdx));
                        end
                    end

                    %smoothed velocity + lick rate
                    if params.binsize_deg <= 2
                        WinSize = 11;
                    else
                        WinSize = 2;
                    end
                    velocCountsSmooth(tr,:) = smooth2005(velocCounts(tr,:),WinSize,'rlowess'); % ASK JLK: why use different smoothing for lick and velocity
                    lickRateSmooth(tr,:) = smoothGaussianMultiple(lickRate(tr,:), 2);

                    %save each trial to the cell after the last zone number
                    % for this zone (e.g., 4th cell in row if 3 reward zones)
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).numRewards = numRewards(tr);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).totalTime = totalTime(tr);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).binEdges = binEdges;
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).lickBehavior = lickBehavior(tr);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).lickCounts = lickCounts(tr,:);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).lickProbability = lickProbability(tr,:);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).lickRate = lickRate(tr,:);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).lickRateSmooth = lickRateSmooth(tr,:);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).velocAvg = velocAvg(tr,:);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).velocCounts = velocCounts(tr,:);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).velocCountsSmooth = velocCountsSmooth(tr,:);
                    statsByRewardTrial{znType,length(znStartDegs)+1}(znNum + (tr-1)*length(znStartDegs)).fileInfo = virmen_fileInfo;

                end%for tr

                lickBehavior = extractStructFields(lickBehavior);

                %save to statsByRewardTrial struct
                statsByRewardTrial{znType,znNum}.numRewards = numRewards;
                statsByRewardTrial{znType,znNum}.totalTime = totalTime;
                statsByRewardTrial{znType,znNum}.binEdges = binEdges;
                statsByRewardTrial{znType,znNum}.lickBehavior = lickBehavior;
                statsByRewardTrial{znType,znNum}.lickCounts = lickCounts;
                statsByRewardTrial{znType,znNum}.lickProbability = lickProbability;
                statsByRewardTrial{znType,znNum}.lickRate = lickRate;
                statsByRewardTrial{znType,znNum}.lickRateSmooth = lickRateSmooth;
                statsByRewardTrial{znType,znNum}.velocAvg = velocAvg;
                statsByRewardTrial{znType,znNum}.velocCounts = velocCounts;
                statsByRewardTrial{znType,znNum}.velocCountsSmooth = velocCountsSmooth;
                statsByRewardTrial{znType,znNum}.fileInfo = virmen_fileInfo;

            else%no data for this zone
                statsByRewardTrial{znType,znNum} = [];
            end%~isempty
        end%for znNum
    end%for znType
end%if length(lapStartIdx) <= 1

%% Get stats across trials throughout session %%
%adapted from getSessionStats_linearJLK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% some initial variables for the trial %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(lapStartIdx) <= 1 %no completed laps

    %save to statsByRewardSession struct
    statsByRewardSession = [];

else %multiple completed laps

    %initialize output struct
    statsByRewardSession = cell(3,1);

    for znType = 1:size(rawDataByTrial,1)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% determine zones %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        if znType == 1%reward
            znStartDegs = params.Rzones;
        elseif znType == 2%control
            znStartDegs = params.NRzones;
        elseif znType == 3%alt control
            znStartDegs = params.AltNRzones;
        end

        statsByRewardSession{znType,1}.binEdges = -params.gapBefore:params.binsize_deg:params.gapAfter;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% reward, time values %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        statsByRewardSession{znType,1}.numRewards = sum([statsByRewardTrial{znType, length(znStartDegs)+1}.numRewards]);
        statsByRewardSession{znType,1}.totalTime = sum([statsByRewardTrial{znType,length(znStartDegs)+1}.totalTime]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% licking behavior %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %initialize variables
        statsByRewardSession{znType,1}.lickBehavior = [];
        tmplicksInZone = [];
        tmplickCounts = [];
        tmplickRateInZone = [];
        tmplickProbability = [];
        tmplickRate = [];

        %get info for each trial
        for tr = 1:length(statsByRewardTrial{znType,length(znStartDegs)+1})
            if isfield(statsByRewardTrial{znType,length(znStartDegs)+1}(tr).lickBehavior, 'licksInZone')
                tmplicksInZone(tr) = statsByRewardTrial{znType,length(znStartDegs)+1}(tr).lickBehavior.licksInZone;
                tmplickCounts(tr,:) = statsByRewardTrial{znType,length(znStartDegs)+1}(tr).lickCounts;
                tmplickRateInZone(tr)  = statsByRewardTrial{znType,length(znStartDegs)+1}(tr).lickBehavior.lickRateInZone;
                tmplickProbability(tr,:) = statsByRewardTrial{znType,length(znStartDegs)+1}(tr).lickProbability;
                tmplickRate(tr,:) = statsByRewardTrial{znType,length(znStartDegs)+1}(tr).lickRate;
            end
        end

        %save info to statsByRewardSession struct
        statsByRewardSession{znType,1}.lickBehavior.licksInZone = sum(tmplicksInZone);
        statsByRewardSession{znType,1}.lickBehavior.licksTotal = sum(tmplickCounts,'all');
        statsByRewardSession{znType,1}.lickBehavior.percentLicksInZone = sum(tmplicksInZone) / sum(tmplickCounts,'all');
        statsByRewardSession{znType,1}.lickBehavior.lickRateInZone = mean(tmplickRateInZone);
        statsByRewardSession{znType,1}.lickBehavior.lickCounts = tmplickCounts;
        statsByRewardSession{znType,1}.lickBehavior.lickCountsSum = sum(tmplickCounts);
        statsByRewardSession{znType,1}.lickBehavior.lickProbability = tmplickProbability;
        statsByRewardSession{znType,1}.lickBehavior.lickProbabilityAvg = mean(tmplickProbability);
        statsByRewardSession{znType,1}.lickBehavior.lickRate = tmplickRate;
        statsByRewardSession{znType,1}.lickBehavior.lickRateAvg = mean(tmplickRate);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% velocity behavior %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %initialize variables
        statsByRewardSession{znType,1}.velocBehavior = [];
        tmpvelocCounts = [];
        tmpvelocCountsSmooth = [];

        %get info for each trial
        for tr = 1:length(statsByRewardTrial{znType,length(znStartDegs)+1})
            if isfield(statsByRewardTrial{znType,length(znStartDegs)+1}(tr).lickBehavior, 'velocCounts')
                tmpvelocCounts(tr,:) = statsByRewardTrial{znType,length(znStartDegs)+1}(tr).velocCounts;
                tmpvelocCountsSmooth(tr,:) = statsByRewardTrial{znType,length(znStartDegs)+1}(tr).velocCountsSmooth;
            end
        end

        %velocity counts total by bin
        statsByRewardSession{znType,1}.velocBehavior.velocCounts = tmpvelocCounts;
        statsByRewardSession{znType,1}.velocBehavior.velocCountsAvg = mean(tmpvelocCounts);
        statsByRewardSession{znType,1}.velocBehavior.velocCountsAvgSmooth = mean(tmpvelocCountsSmooth);

    end%znType

end%if length(lapStartIdx) <= 1

%% Save stats %%
filename = [saveBehaviorPath '\' 'statsByRewardTrial.mat'];
filename2 = [saveBehaviorPath '\' 'statsByRewardSession.mat'];
save(filename, 'statsByRewardTrial');
save(filename2, 'statsByRewardSession');
if length(lapStartIdx) > 1
    filename3 = [saveBehaviorPath '\' 'rawDataByTrial.mat'];
    save(filename3, 'rawDataByTrial');
end
disp(['saved trial stats for ' saveBehaviorPath])
end%function