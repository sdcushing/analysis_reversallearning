function [rawDataByLap, statsByLap, statsBySession] = getLapBehaviorStats_linearJLK(rawDataBySession, virmen_fileInfo, params, saveBehaviorPath)
%Note: only completed laps

%% Get stats by lap %%
%adapted from getTrialByTrialStats_linearJLK

% split data up into individual laps (0 to 360 deg)
%lapStartIdx = index of starting a new lap, going from 0 to 360 deg
lapStartIdx = find(diff(rawDataBySession.currentDeg)<-300)+1;
lapStartIdx = [1;lapStartIdx(1:end)];

% Calculate lap by lap stats
clear numRewards totalTime lickBehavior lickCounts lickProbability lickRate velocAvg velocCounts velocCountsSmooth rawDataByLap statsByLap statsBySession

%split data up into individual laps based on completion of track
if length(lapStartIdx) <= 1

    %save to rawDataByLap and statsByLap structs
    rawDataByLap = [];
    statsByLap = [];

else

    for ii = 1:size(lapStartIdx,1)-1
        rawDataByLap(ii) = structfun(@(x) x(lapStartIdx(ii):lapStartIdx(ii+1)-1), rawDataBySession, 'UniformOutput', false);
    end

    %initialize variables
    histEdges = 0:params.binsize_deg:360;
    lickCounts = zeros(length(rawDataByLap),length(histEdges)-1);
    lickProbability = zeros(length(rawDataByLap),length(histEdges)-1);
    lickRate = zeros(length(rawDataByLap),length(histEdges)-1);
    velocCounts = zeros(length(rawDataByLap),length(histEdges)-1);
    velocCountsSmooth = zeros(length(rawDataByLap),length(histEdges)-1);

    for ii = 1:length(rawDataByLap)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% some initial variables for the lap %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lickTimes = [0; diff(rawDataByLap(ii).licks)];

        timeout = find(rawDataByLap(ii).currentZone == params.timeoutZone);
        %find times during timeout periods to exclude
        if ~isempty(timeout)
            lickTimes(timeout) = 0; %converts any licks that occurred during timeout to zero
        end
        totalLicks = sum(lickTimes); %total lickcount per lap

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% reward, time values %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numRewards(ii,:) = max(rawDataByLap(ii).rewards)-min(rawDataByLap(ii).rewards);
        totalTime(ii,:) = rawDataByLap(ii).vrTime(end)-rawDataByLap(ii).vrTime(1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% licking behavior %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %adapted from getLickInfoForTrial_linearJLK

        %initialize variables
        licksInReward = zeros(1, length(params.Rzones));
        licksInAnticip = zeros(1, length(params.Rzones));
        licksInControl = zeros(1, length(params.NRzones));
        licksInNevControl = zeros(1, length(params.NevRzones));

        timeInReward = zeros(1, length(params.Rzones));
        timeInAnticip = zeros(1, length(params.Rzones));
        timeInControl = zeros(1, length(params.NRzones));
        timeInNevControl = zeros(1, length(params.NevRzones));


        %number of licks in anticipatory and reward zone
        for rz = 1:length(params.Rzones)
            licksInReward(rz) = sum(lickTimes(find(rawDataByLap(ii).currentDeg >= params.Rzones(rz) & rawDataByLap(ii).currentDeg < wrapTo360(params.Rzones(rz) + virmen_fileInfo.cueSize))));
            licksInAnticip(rz) = sum(lickTimes(find(rawDataByLap(ii).currentDeg >= wrapTo360(params.Rzones(rz) - virmen_fileInfo.cueSize ) & rawDataByLap(ii).currentDeg < params.Rzones(rz) )));

            timeInReward(rz) = sum(diff(rawDataByLap(ii).vrTime(find(rawDataByLap(ii).currentDeg >= params.Rzones(rz) & rawDataByLap(ii).currentDeg < wrapTo360(params.Rzones(rz) + virmen_fileInfo.cueSize )))));
            timeInAnticip(rz) = sum(diff(rawDataByLap(ii).vrTime(find(rawDataByLap(ii).currentDeg >= wrapTo360(params.Rzones(rz) - virmen_fileInfo.cueSize ) & rawDataByLap(ii).currentDeg < params.Rzones(rz) ))));
        end

        %number of licks in control zone (could be different number than reward zones)
        for nrz = 1:length(params.NRzones)
            licksInControl(nrz) = sum(lickTimes(find(rawDataByLap(ii).currentDeg >= params.NRzones(nrz) & rawDataByLap(ii).currentDeg < wrapTo360(params.NRzones(nrz) + virmen_fileInfo.cueSize))));
            timeInControl(nrz) = sum(diff(rawDataByLap(ii).vrTime(find(rawDataByLap(ii).currentDeg >= params.NRzones(nrz) & rawDataByLap(ii).currentDeg < wrapTo360(params.NRzones(nrz) + virmen_fileInfo.cueSize)))));
        end

        %number of licks in alt control zone (could be different number than reward zones)
        if ~isnan(params.NevRzones(1))
            for altnrz = 1:length(params.NevRzones)
                licksInNevControl(altnrz) = sum(lickTimes(find(rawDataByLap(ii).currentDeg >= params.NevRzones(altnrz) & rawDataByLap(ii).currentDeg < wrapTo360(params.NevRzones(altnrz) + virmen_fileInfo.cueSize))));
                timeInNevControl(altnrz) = sum(diff(rawDataByLap(ii).vrTime(find(rawDataByLap(ii).currentDeg >= params.NevRzones(altnrz) & rawDataByLap(ii).currentDeg < wrapTo360(params.NevRzones(altnrz) + virmen_fileInfo.cueSize)))));
            end
        end

        %number of licks per zone
        lickBehavior(ii).licksInReward = sum(licksInReward);
        lickBehavior(ii).licksInAnticip = sum(licksInAnticip);
        lickBehavior(ii).licksInControl = sum(licksInControl);
        lickBehavior(ii).licksInNevControl = sum(licksInNevControl);


        %number of licks per zone
        lickBehavior(ii).timeInReward = sum(timeInReward);
        lickBehavior(ii).timeInAnticip = sum(timeInAnticip);
        lickBehavior(ii).timeInControl = sum(timeInControl);
        lickBehavior(ii).timeInNevControl = sum(timeInNevControl);


        %lick rate per zone
        lickBehavior(ii).lickRateInReward = lickBehavior(ii).licksInReward/lickBehavior(ii).timeInReward;
        lickBehavior(ii).lickRateInAnticip = lickBehavior(ii).licksInAnticip/lickBehavior(ii).timeInAnticip;
        lickBehavior(ii).lickRateInControl = lickBehavior(ii).licksInControl/lickBehavior(ii).timeInControl;
        lickBehavior(ii).lickRateInNevControl = lickBehavior(ii).licksInNevControl/lickBehavior(ii).timeInNevControl;


        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% licking counts %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %adapted from lickCountsForTrial_linearNJ

        %calculate lick bins across the lap
        for i = 1:length(histEdges)-1
            binIdx = find(rawDataByLap(ii).currentDeg >= histEdges(i) & rawDataByLap(ii).currentDeg < histEdges(i+1));
            if ~isempty(binIdx)
                lickCounts(ii,i) = sum(lickTimes(binIdx));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% licking probability %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %adapted from lickProbabilityForTrial_linearJLK

        %calculate lick bins across the lap
        for i = 1:length(histEdges)-1
            binIdx = find(rawDataByLap(ii).currentDeg >= histEdges(i) & rawDataByLap(ii).currentDeg < histEdges(i+1));
            if ~isempty(binIdx)
                lickProbability(ii,i) = sum(lickTimes(binIdx))/totalLicks;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%
        %%%%% lick rate %%%%%
        %%%%%%%%%%%%%%%%%%%%%
        %adapted from lickRateForTrial_linearJLK

        %calculate lick bins across the lap
        for i = 1:length(histEdges)-1
            binIdx = find(rawDataByLap(ii).currentDeg >= histEdges(i) & rawDataByLap(ii).currentDeg < histEdges(i+1));
            if ~isempty(binIdx)
                binTime = rawDataByLap(ii).vrTime(binIdx);
                lickRate(ii,i) = sum(lickTimes(binIdx)) / (binTime(end)-binTime(1));
            else
                lickRate(ii,i) = 0;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% velocity average %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        velocAvg(ii,:) = (wrapTo360(rawDataByLap(ii).currentDeg(end)-rawDataByLap(ii).currentDeg(1)))/totalTime(ii);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% velocity counts %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        WinSize = 11;

        %calculate instantaneous velocity from smoothed coordinates across lap
        velocVect = zeros(1,length(rawDataByLap(ii).currentDeg));

        for i = 1:length(rawDataByLap(ii).currentDeg)-1
            totaldist = abs(rawDataByLap(ii).currentDeg(i+1) - rawDataByLap(ii).currentDeg(i));
            totaltime = (rawDataByLap(ii).vrTime(i+1) - rawDataByLap(ii).vrTime(i));
            if rawDataByLap(ii).currentZone(i+1) == params.timeoutZone
                % exclude time-out periods from velocity calculation
                velocVect(i) = 0; %timeout velocity printed as zero (no actual movement through VR during timeout period)
            else
                velocVect(i) = totaldist / totaltime;
            end
        end

        %get velocity bins
        for i = 1:length(histEdges)-1
            binIdx = find(rawDataByLap(ii).currentDeg >= histEdges(i) & rawDataByLap(ii).currentDeg < histEdges(i+1));
            if ~isempty(binIdx)
                velocCounts(ii,i) = mean(velocVect(binIdx));
            end
        end

        %smoothed velocity
        velocCountsSmooth(ii,:) = smooth2005(velocCounts(ii,:),WinSize,'rlowess');
    end

    lickRateSmooth = smoothGaussianMultiple(lickRate, 2);

    %save to statsByLap struct
    statsByLap.numRewards = numRewards;
    statsByLap.totalTime = totalTime;
    statsByLap.binEdges = 0:params.binsize_deg:360;
    statsByLap.lickBehavior = lickBehavior;
    statsByLap.lickCounts = lickCounts;
    statsByLap.lickProbability = lickProbability;
    statsByLap.lickRate = lickRate;
    statsByLap.lickRateSmooth = lickRateSmooth;
    statsByLap.velocAvg = velocAvg;
    statsByLap.velocCounts = velocCounts;
    statsByLap.velocCountsSmooth = velocCountsSmooth;
    statsByLap.fileInfo = virmen_fileInfo;
end%if length(lapStartIdx) <= 1

%% Get stats across laps throughout session %%
%adapted from getSessionStats_linearJLK

if ~isempty(statsByLap)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% some initial variables for the trial %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    statsBySession.numRewards = sum(statsByLap.numRewards);
    statsBySession.totalTime = sum(statsByLap.totalTime);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% licking behavior %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %sum lap info for the whole session
    statsBySession.lickBehavior.licksInReward = sum([statsByLap.lickBehavior.licksInReward]);
    statsBySession.lickBehavior.licksInAnticip = sum([statsByLap.lickBehavior.licksInAnticip]);
    statsBySession.lickBehavior.licksInControl = sum([statsByLap.lickBehavior.licksInControl]);
    statsBySession.lickBehavior.licksInNevControl = sum([statsByLap.lickBehavior.licksInNevControl]);

    %calculate percentages
    statsBySession.lickBehavior.licksTotal = sum(statsByLap.lickCounts, 'all');
    statsBySession.lickBehavior.percentLicksInReward = (statsBySession.lickBehavior.licksInReward/statsBySession.lickBehavior.licksTotal)*100;
    statsBySession.lickBehavior.percentLicksInAnticip = (statsBySession.lickBehavior.licksInAnticip/statsBySession.lickBehavior.licksTotal)*100;
    statsBySession.lickBehavior.percentLicksInControl = (statsBySession.lickBehavior.licksInControl/statsBySession.lickBehavior.licksTotal)*100;
    statsBySession.lickBehavior.percentLicksInNevControl = (statsBySession.lickBehavior.licksInNevControl/statsBySession.lickBehavior.licksTotal)*100;

    %lick count total by bin
    statsBySession.lickBehavior.lickCountsSum = sum(statsByLap.lickCounts);

    %lick probability by bin
    statsBySession.lickBehavior.lickProbabilityAvg = mean(statsByLap.lickProbability);

    %lick rate by bin
    statsBySession.lickBehavior.lickRateAvg = mean(statsByLap.lickRate);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% velocity behavior %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %velocity counts total by bin
    statsBySession.velocBehavior.velocCountsSum =  sum(statsByLap.velocCounts);
    statsBySession.velocBehavior.velocCountsSmoothSum =  sum(statsByLap.velocCountsSmooth);

    %velocity counts average by bin
    statsBySession.velocBehavior.velocCountsAvg = mean(statsByLap.velocCounts);
    statsBySession.velocBehavior.velocCountsSmoothAvg = mean(statsByLap.velocCountsSmooth);


else%if not one full lap
    statsBySession = [];
end

%% Save data and stats %%
filename = [saveBehaviorPath '\' 'fileInfo.mat'];
filename2 = [saveBehaviorPath '\' 'statsByLap.mat'];
filename3 = [saveBehaviorPath '\' 'statsBySession.mat'];
filename4 = [saveBehaviorPath '\' 'rawDataByLap.mat'];
save(filename, 'virmen_fileInfo');
save(filename2, 'statsByLap');
save(filename3, 'statsBySession');
save(filename4, 'rawDataByLap');

end%function