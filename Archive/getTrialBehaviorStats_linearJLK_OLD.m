function [rawDataByTrial, statsByRewardTrial, statsByRewardSession] = getTrialBehaviorStats_linearJLK(virmen_dataStruct, virmen_fileInfo, params, saveBehaviorPath)
%Note: only completed laps

%% GET STATS BY TRIAL
%adapted from getTrialByTrialStats_linearJLK

% split data up into individual laps (0 to 360 deg)
%lapStartIdx = index of starting a new lap, going from 0 to 360 deg
lapStartIdx = find(diff(virmen_dataStruct.currentDeg)<-300)+1;
lapStartIdx = [1;lapStartIdx(1:end)];

% Calculate trial by trial stats
clear numRewards totalTime lickBehavior lickCounts lickProbability lickRate velocAvg velocCounts velocCountsSmooth statsByRewardTrial statsByRewardSession

if length(lapStartIdx) < 2 %no completed laps
            
    %save to statsByRewardTrial struct
    statsByRewardTrial = [];

else %at least one completed lap
   
    %parameters
    endBin = abs(params.gapBefore) + abs(params.gapAfter);
    binEdges = -params.gapBefore:params.binsize_deg:params.gapAfter;
    histEdges = 0:params.binsize_deg:endBin;
    
    %padded zones are X degrees (as defined by gap) before and after each reward zone
    paddedZones(:,1) = wrapTo360(params.Rzones - params.gapBefore);
    paddedZones(:,2) = wrapTo360(params.Rzones + virmen_fileInfo.cueSize + params.gapAfter);
    
    %split data up into individual trials based on reward zones
    allC = [];
    for p = 1:length(params.Rzones) %in the order of reward zones separately (1,2,3,4)
        temp = isExcluded(virmen_dataStruct.currentDeg, paddedZones(p,:));
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
                idx2excl = find(wrapTo360(virmen_dataStruct.currentDeg(C(:,2)) - virmen_dataStruct.currentDeg(C(:,1))) < params.gapBefore + params.gapAfter);
                C(idx2excl,:) = [];
                allC = [allC; C];
            else
                C = [PZ(1:2:end-1),PZ(2:2:end-1)];
                idx2excl = find(wrapTo360(virmen_dataStruct.currentDeg(C(:,2)) - virmen_dataStruct.currentDeg(C(:,1))) < params.gapBefore + params.gapAfter);
                if ~isempty(idx2excl)
                    C(idx2excl,:) = [];
                end
                allC = [allC; C];
            end
            for ii = 1:size(C,1)
                rawDataByTrial{p}(ii) = structfun(@(x) x(C(ii,1):C(ii,2)), virmen_dataStruct, 'UniformOutput', false);
            end
            clear PZ
        end
        allC = sort(allC,1); %last struct uses all 4 reward zones in order of appaearance
        
        for ii = 1:size(allC,1)
            rawDataByTrial{length(paddedZones)+1}(ii) = structfun(@(x) x(allC(ii,1):allC(ii,2)), virmen_dataStruct, 'UniformOutput', false);
        end
        
    end
    
    % stats for reward trials
    for t = 1:length(rawDataByTrial)
        if length(rawDataByTrial{t}) > 0 %added this 3/4/22 to deal with rare cases where there was no reward trial

            %initialize variables
            lickCounts = zeros(length(rawDataByTrial{t}),length(histEdges)-1);
            lickProbability = zeros(length(rawDataByTrial{t}),length(histEdges)-1);
            lickRate = zeros(length(rawDataByTrial{t}),length(histEdges)-1);
            velocCounts = zeros(length(rawDataByTrial{t}),length(histEdges)-1);
            velocCountsSmooth = zeros(length(rawDataByTrial{t}),length(histEdges)-1);

            for ii = 1:length(rawDataByTrial{t})

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% some initial variables for the trial %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                lickTimes = [0; diff(rawDataByTrial{t}(ii).licks)];

                timeout = find(rawDataByTrial{t}(ii).currentZone == params.timeoutZone);
                %find times during timeout periods to exclude
                if ~isempty(timeout)
                    lickTimes(timeout) = 0; %converts any licks that occurred during timeout to zero
                end
                totalLicks = sum(lickTimes); %total lickcount per trial

                vel = rawDataByTrial{t}(ii).currentDeg;
                vel = [0; wrapTo360(diff(vel))];
                vel(vel>300) = 0;
                cumvel = cumsum(vel); %final vel = degrees that go from zero to 70

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% reward, time values %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                numRewards(ii,:) = max(rawDataByTrial{t}(ii).rewards)-min(rawDataByTrial{t}(ii).rewards);
                totalTime(ii,:) = rawDataByTrial{t}(ii).vrTime(end)-rawDataByTrial{t}(ii).vrTime(1);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% licking behavior %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %adapted from getLickInfoForRewardTrial_linearJLK

                %number of licks in each zone
                currRZ = lookup2(rawDataByTrial{t}(ii).currentDeg(1), wrapTo360(params.Rzones - params.gapBefore));
                currNRZ = lookup2(rawDataByTrial{t}(ii).currentDeg(1), wrapTo360(params.NRzones - params.gapBefore));
                licksInReward = sum(lickTimes(find(rawDataByTrial{t}(ii).currentDeg >= params.Rzones(currRZ) & rawDataByTrial{t}(ii).currentDeg < wrapTo360(params.Rzones(currRZ) + virmen_fileInfo.cueSize))));
                licksInAnticip = sum(lickTimes(find(rawDataByTrial{t}(ii).currentDeg >= wrapTo360(params.Rzones(currRZ) - virmen_fileInfo.cueSize ) & rawDataByTrial{t}(ii).currentDeg < params.Rzones(currRZ) )));
                licksInControl = sum(lickTimes(find(rawDataByTrial{t}(ii).currentDeg >= params.NRzones(currNRZ) & rawDataByTrial{t}(ii).currentDeg < wrapTo360(params.NRzones(currNRZ) + virmen_fileInfo.cueSize))));

                timeInReward = sum(diff(rawDataByTrial{t}(ii).vrTime(find(rawDataByTrial{t}(ii).currentDeg >= params.Rzones(currRZ) & rawDataByTrial{t}(ii).currentDeg < wrapTo360(params.Rzones(currRZ) + virmen_fileInfo.cueSize )))));
                timeInAnticip = sum(diff(rawDataByTrial{t}(ii).vrTime(find(rawDataByTrial{t}(ii).currentDeg >= wrapTo360(params.Rzones(currRZ) - virmen_fileInfo.cueSize ) & rawDataByTrial{t}(ii).currentDeg < params.Rzones(currRZ) ))));
                timeInControl = sum(diff(rawDataByTrial{t}(ii).vrTime(find(rawDataByTrial{t}(ii).currentDeg >= params.NRzones(currNRZ) & rawDataByTrial{t}(ii).currentDeg < wrapTo360(params.NRzones(currNRZ) + virmen_fileInfo.cueSize)))));

                %lick rate per zone
                lickRateInReward = licksInReward/timeInReward;
                lickRateInAnticip = licksInAnticip/timeInAnticip;
                lickRateInControl = licksInControl/timeInControl;

                %get lick outs
                lickBehavior(ii).licksInReward = licksInReward;
                lickBehavior(ii).licksInAnticip = licksInAnticip;
                lickBehavior(ii).licksInControl = licksInControl;

                lickBehavior(ii).timeInReward = timeInReward;
                lickBehavior(ii).timeInAnticip = timeInAnticip;
                lickBehavior(ii).timeInControl = timeInControl;

                lickBehavior(ii).lickRateInReward = lickRateInReward;
                lickBehavior(ii).lickRateInAnticip = lickRateInAnticip;
                lickBehavior(ii).lickRateInControl = lickRateInControl;

                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% licking counts %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                %adapted from lickCountsForRewardTrial_linearJLK

                %calculate lick bins across the trial
                for i = 1:length(histEdges)-1
                    binIdx = find(cumvel >= histEdges(i) & cumvel < histEdges(i+1));
                    if ~isempty(binIdx)
                        lickCounts(ii,i) = sum(lickTimes(binIdx));
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% licking probability %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %adapted from lickProbabilityForRewardTrial_linearJLK

                %calculate lick bins across the trial
                for i = 1:length(histEdges)-1
                    binIdx = find(cumvel >= histEdges(i) & cumvel < histEdges(i+1));
                    if ~isempty(binIdx)
                        lickProbability(ii,i) = sum(lickTimes(binIdx))/totalLicks;
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%
                %%%%% lick rate %%%%%
                %%%%%%%%%%%%%%%%%%%%%
                %adapted from lickRateForTrial_linearJLK

                %calculate lick bins across the trial
                for i = 1:length(histEdges)-1
                    binIdx = find(cumvel >= histEdges(i) & cumvel < histEdges(i+1));
                    if ~isempty(binIdx)
                        binTime = rawDataByTrial{t}(ii).vrTime(binIdx);
                        lickRate(ii,i) = sum(lickTimes(binIdx)) / (binTime(end)-binTime(1));
                    else
                        lickRate(ii,i) = 0;
                    end
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% velocity average %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                velocAvg(ii,:) = sum(wrapTo360(diff(rawDataByTrial{t}(ii).currentDeg)))/totalTime(ii);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% velocity counts %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                WinSize = 11;

                %calculate instantaneous velocity from smoothed coordinates across trial
                velocVect = zeros(1,length(rawDataByTrial{t}(ii).currentDeg));
                timeDiff = [0; diff(rawDataByTrial{t}(ii).vrTime)];

                for i = 1:length(rawDataByTrial{t}(ii).currentDeg)
                    if rawDataByTrial{t}(ii).currentZone(i) == params.timeoutZone
                        % exclude time-out periods from velocity calculation
                        velocVect(i) = 0; %timeout velocity printed as zero (no actual movement through VR during timeout period)
                    else
                        velocVect(i) = vel(i) ./ timeDiff(i);
                    end
                end

                %get velocity bins
                for i = 1:length(histEdges)-1
                    binIdx = find(cumvel >= histEdges(i) & cumvel < histEdges(i+1));
                    if ~isempty(binIdx)
                        velocCounts(ii,i) = mean(velocVect(binIdx));
                    end
                end

                %smoothed velocity
                velocCountsSmooth(ii,:) = smooth2005(velocCounts(ii,:),WinSize,'rlowess');

            end%for ii

            lickRateSmooth = smoothGaussianMultiple(lickRate, 2);
            lickBehavior = extractStructFields(lickBehavior);
            
            %save to statsByRewardTrial struct
            statsByRewardTrial{t}.numRewards = numRewards;
            statsByRewardTrial{t}.totalTime = totalTime;
            statsByRewardTrial{t}.binEdges = binEdges;
            statsByRewardTrial{t}.lickBehavior = lickBehavior;
            statsByRewardTrial{t}.lickCounts = lickCounts;
            statsByRewardTrial{t}.lickProbability = lickProbability;
            statsByRewardTrial{t}.lickRate = lickRate;
            statsByRewardTrial{t}.lickRateSmooth = lickRateSmooth;
            statsByRewardTrial{t}.velocAvg = velocAvg;
            statsByRewardTrial{t}.velocCounts = velocCounts;
            statsByRewardTrial{t}.velocCountsSmooth = velocCountsSmooth;
            statsByRewardTrial{t}.fileInfo = virmen_fileInfo;

        end%if length(rawDataByTrial{t}) > 0
    end%for t
end%if length(lapStartIdx) <= 1

%% GET STATS BY ACROSS TRIALS THROUGHOUT SESSION
%adapted from getSessionStats_linearJLK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% some initial variables for the trial %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(lapStartIdx) <= 1 %no completed laps

    %save to statsByRewardSession struct
    statsByRewardSession = [];

else %multiple completed laps
     
    statsByRewardSession.binEdges = -params.gapBefore:params.binsize_deg:params.gapAfter;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% reward, time values %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    statsByRewardSession.numRewards = sum(statsByRewardTrial{end}.numRewards);
    statsByRewardSession.totalTime = sum(statsByRewardTrial{end}.totalTime);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% licking behavior %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %sum trial info for the whole session
    statsByRewardSession.lickBehavior.licksInReward = sum([statsByRewardTrial{end}.lickBehavior.licksInReward]);
    statsByRewardSession.lickBehavior.licksInAnticip = sum([statsByRewardTrial{end}.lickBehavior.licksInAnticip]);
    statsByRewardSession.lickBehavior.licksInControl = sum([statsByRewardTrial{end}.lickBehavior.licksInControl]);
    
    %calculate percentages
    statsByRewardSession.lickBehavior.licksTotal = sum(statsByRewardTrial{end}.lickCounts, 'all');
    statsByRewardSession.lickBehavior.percentLicksInReward = (statsByRewardSession.lickBehavior.licksInReward/statsByRewardSession.lickBehavior.licksTotal)*100;
    statsByRewardSession.lickBehavior.percentLicksInAnticip = (statsByRewardSession.lickBehavior.licksInAnticip/statsByRewardSession.lickBehavior.licksTotal)*100;
    statsByRewardSession.lickBehavior.percentLicksInControl = (statsByRewardSession.lickBehavior.licksInControl/statsByRewardSession.lickBehavior.licksTotal)*100;

    %calculate rates
    statsByRewardSession.lickBehavior.lickRateInReward = mean([statsByRewardTrial{end}.lickBehavior.licksInReward]);
    statsByRewardSession.lickBehavior.lickRateInAnticip = mean([statsByRewardTrial{end}.lickBehavior.licksInAnticip]);
    statsByRewardSession.lickBehavior.lickRateInControl = mean([statsByRewardTrial{end}.lickBehavior.licksInControl]);

    %lick counts total by bin
    statsByRewardSession.lickBehavior.lickCountsSum = sum([statsByRewardTrial{end}.lickCounts]);
    
    %lick probability by bin
    statsByRewardSession.lickBehavior.lickProbabilityAvg = mean([statsByRewardTrial{end}.lickProbability]);

    %lick rate
    statsByRewardSession.lickBehavior.lickRateAvg = mean([statsByRewardTrial{end}.lickRate]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% velocity behavior %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %velocity counts total by bin
    statsByRewardSession.velocBehavior.velocCountsAvg = mean([statsByRewardTrial{end}.velocCounts]);
    statsByRewardSession.velocBehavior.velocCountsAvgSmooth = mean([statsByRewardTrial{end}.velocCountsSmooth]);

end%if length(lapStartIdx) <= 1

%% Save data and stats
filename = [saveBehaviorPath '\' 'statsByRewardTrial'];
filename2 = [saveBehaviorPath '\' 'statsByRewardSession'];
if ~exist(filename) || params.rewrite.behavior
    save(filename, 'statsByRewardTrial');
    save(filename2, 'statsByRewardSession');
end

end%function