function [rawDataBySessionNeural] = getRipplesTmp_DC(dirs, params, saveNeuralPath, plotRipples)
%adapted from filtereeg2_Intan.m, extractripples3.m,
% ripplepostfileprocess2, findpowerratioripplevsabove2
%last check JLK 10/29/25

%% extract ripples across session %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load session data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([saveNeuralPath '\rawDataBySessionNeural.mat'])
load([saveNeuralPath '\sessionPyrLayerInfo.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% initialize params and load filters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = [];
beta = [];
delta = [];
tdbratio = [];
ripple = [];
ripples = [];
ripplesBad = [];
chCtr = 0;
layerChans = sessionPyrLayerInfo.pyrLayerCA1;%CA1 only
%layerChans = [sessionPyrLayerInfo.pyrLayerCA3 sessionPyrLayerInfo.pyrLayerCA1];%CA1+CA3
load([dirs.code 'ripplefilter.mat'])
load([dirs.code 'thetafilter.mat'])
load([dirs.code 'betafilter.mat'])
load([dirs.code 'deltafilter.mat'])

%%%%% add filter description and parameters %%%%%
ripple.descript = ripplefilter.descript;
ripple.kernel = ripplefilter.kernel;
ripple.samprate = ripplefilter.samprate;
smoothing_width_rip = 0.004; % 4 ms
smoothing_kernel_rip = gaussian(smoothing_width_rip*ripple.samprate, ceil(8*smoothing_width_rip*ripple.samprate));
smoothing_width_tdb = 1; % 1 s
smoothing_kernel_tdb = gaussian(smoothing_width_tdb*ripple.samprate, ceil(8*smoothing_width_tdb*ripple.samprate));
minRipDur = round(params.ripple.minRipDur * ripple.samprate);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find ripples per channel %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ch = layerChans
    chCtr = chCtr+1;
    %ripple filter (150-250 Hz)
    ripple.filtdata(chCtr,:) = filtfilt(ripplefilter.kernel, 1 , rawDataBySessionNeural.lfpData(ch,:));

    %Hilbert transform ripple data
    hdata = hilbert(ripple.filtdata(chCtr,:));
    ripple.phase(chCtr,:) = angle(hdata);
    ripple.env(chCtr,:) = abs(hdata);
    ripple.env(chCtr,:) = smoothvect(ripple.env(chCtr,:), smoothing_kernel_rip);%smooth env
    clear hdata

    for t = 1:length(params.ripple.nstdEnv)%loop through detection thresholds

        %find baseline, std, and thresh for this channel's ripple env
        envBaseline = mean(ripple.env(chCtr,:));
        envStd = std(ripple.env(chCtr,:));
        envThresh = envBaseline + params.ripple.nstdEnv(t) * envStd;

        %extract the events if this is a valid trace
        if (envThresh > 0) && any(find(ripple.env(chCtr,:)<envBaseline))
            %find ripple events
            tmprip = [];
            tmprip = extractevents(ripple.env(chCtr,:), envThresh, envBaseline, 0, minRipDur, 0)';
            %start, middle (of energy), and end indices and times
            ripples(chCtr,t).startind = tmprip(:,1);
            ripples(chCtr,t).midind = tmprip(:,8);
            ripples(chCtr,t).endind = tmprip(:,2);
            ripples(chCtr,t).starttime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).startind / ripple.samprate;
            ripples(chCtr,t).midtime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).midind / ripple.samprate;
            ripples(chCtr,t).endtime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).endind / ripple.samprate;
            %ripple characteristics
            ripples(chCtr,t).peak = tmprip(:,3);
            ripples(chCtr,t).energy = tmprip(:,7);
            ripples(chCtr,t).maxthresh = (tmprip(:,9) - envBaseline) / envStd;
        else
            ripples(chCtr,t).startind = [];
            ripples(chCtr,t).midind = [];
            ripples(chCtr,t).endind = [];
            ripples(chCtr,t).starttime = [];
            ripples(chCtr,t).midtime = [];
            ripples(chCtr,t).endtime = [];
            ripples(chCtr,t).peak = [];
            ripples(chCtr,t).energy = [];
            ripples(chCtr,t).maxthresh = [];
        end%> thresh + baseline

        %add other info to ripples struct
        ripples(chCtr,t).timeRange = [0 length(ripple.env(chCtr,:))/ripple.samprate] + rawDataBySessionNeural.lfpTime(1);
        ripples(chCtr,t).samprate = ripple.samprate;
        ripples(chCtr,t).envThreshold = envThresh;
        ripples(chCtr,t).envBaseline = envBaseline;
        ripples(chCtr,t).envStd = envStd;
        ripples(chCtr,t).minimum_duration = params.ripple.minRipDur;
        ripples(chCtr,t).chan = ch;

    end%thresholds

end%ch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% collect ripples across channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% concatenate ripples across channels %%%%%
allRipMidInd3SD = []; allRipMidInd5SD = []; allRipCh3SD = []; allRipCh5SD = [];
for ch = 1:length(layerChans)
    %3 SD Threshold
    allRipMidInd3SD = [allRipMidInd3SD; ripples(ch,1).midind];
    allRipCh3SD = [allRipCh3SD; transpose(ones(1,length(ripples(ch,1).midind)))*ch];
    %5 SD Threshold
    allRipMidInd5SD = [allRipMidInd5SD; ripples(ch,2).midind];
    allRipCh5SD = [allRipCh5SD; transpose(ones(1,length(ripples(ch,2).midind)))*ch];
end%ch

%%%%% sort the ripples %%%%%
[allRipMidInd3SD, I3SD] = sort(allRipMidInd3SD);
allRipCh3SD = allRipCh3SD(I3SD);
[allRipMidInd5SD, I5SD] = sort(allRipMidInd5SD);
allRipCh5SD = allRipCh5SD(I5SD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find good ripples across channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% reduce to unique ripples detected at 5 SD %%%%%
allUniqRipMidInd5SD = allRipMidInd5SD(diff(allRipMidInd5SD)>(params.ripple.minRipDur * params.lfp_samprate_down * 4));
allUniqRipCh5SD = allRipCh5SD(diff(allRipMidInd5SD)>(params.ripple.minRipDur * params.lfp_samprate_down * 4));

%%%%% find ripples detected at 5 SD close to ripples detected at 3 SD on 2+ other channels %%%%%
allUniqRipMidInd5SDGood = []; allUniqRipCh5SDGood = []; goodRipCtr = 0;
for r = 1:length(allUniqRipMidInd5SD)
    if sum(allRipMidInd3SD<allUniqRipMidInd5SD(r)+(params.ripple.minRipDur * params.lfp_samprate_down * 4) &...
            allRipMidInd3SD>allUniqRipMidInd5SD(r)-(params.ripple.minRipDur * params.lfp_samprate_down * 4)) >= 3%3 because this ripple + 2 more
        goodRipCtr = goodRipCtr + 1;
        allUniqRipMidInd5SDGood(goodRipCtr) = allUniqRipMidInd5SD(r);
        allUniqRipCh5SDGood(goodRipCtr) = allUniqRipCh5SD(r);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% collect all info for good ripples %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ripplesGood = [];
goodRipCtrPerCh = zeros(1,length(layerChans));
for r = 1:length(allUniqRipMidInd5SDGood)
    goodRipCtrPerCh(allUniqRipCh5SDGood(r)) = goodRipCtrPerCh(allUniqRipCh5SDGood(r))+1;
    goodRip = find(ripples(allUniqRipCh5SDGood(r),2).midind==allUniqRipMidInd5SDGood(r));
    ripplesGood(allUniqRipCh5SDGood(r)).startind(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).startind(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).midind(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).midind(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).endind(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).endind(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).starttime(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).starttime(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).midtime(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).midtime(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).endtime(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).endtime(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).peak(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).peak(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).energy(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).energy(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).maxthresh(goodRipCtrPerCh(allUniqRipCh5SDGood(r))) = ripples(allUniqRipCh5SDGood(r),2).maxthresh(goodRip);
    ripplesGood(allUniqRipCh5SDGood(r)).timeRange = ripples(allUniqRipCh5SDGood(r),2).timeRange;
    ripplesGood(allUniqRipCh5SDGood(r)).samprate = ripples(allUniqRipCh5SDGood(r),2).samprate;
    ripplesGood(allUniqRipCh5SDGood(r)).envThreshold = ripples(allUniqRipCh5SDGood(r),2).envThreshold;
    ripplesGood(allUniqRipCh5SDGood(r)).envBaseline = ripples(allUniqRipCh5SDGood(r),2).envBaseline;
    ripplesGood(allUniqRipCh5SDGood(r)).envStd = ripples(allUniqRipCh5SDGood(r),2).envStd;
    ripplesGood(allUniqRipCh5SDGood(r)).minimum_duration = ripples(allUniqRipCh5SDGood(r),2).minimum_duration;
    ripplesGood(allUniqRipCh5SDGood(r)).chan = ripples(allUniqRipCh5SDGood(r),2).chan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% appply ripple criteria (optional) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ripplesBad = [];
if params.ripple.applyCriteria
    chCtr = 0;
    for ch = layerChans
        chCtr = chCtr+1;
        if any(ripplesGood(chCtr).midind)%if channel has ripples

            %find time around ripple middle (2 s total)
            ripperiods = [ripplesGood(chCtr).midind'-params.ripple.timeAroundRip*ripple.samprate, ripplesGood(chCtr).midind'+params.ripple.timeAroundRip*ripple.samprate];

            %initilize variables
            highFreqVals = [];
            highFreqRatio = [];
            meantdb = [];
            excludeForHighFreq = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForMeanTDB = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForSpeed = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForRawDataNoise = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForNoSpikesPresent = zeros(1,length(ripplesGood(chCtr).midind));
            ripSpikes = nan(100,length(rawDataBySessionNeural.apData),50);

            %find theta/delta+beta ratio for this channel
            if params.ripple.applyTDBRatio
                % theta filter (6-10 Hz)
                theta.filtdata(chCtr,:) = filtfilt(thetafilter.tf.num, 1 , rawDataBySessionNeural.lfpData(ch,:));
                % beta filter (12-30 Hz)
                beta.filtdata(chCtr,:) = filtfilt(betafilter.tf.num, 1 , rawDataBySessionNeural.lfpData(ch,:));
                % delta filter (1-4 Hz)
                delta.filtdata(chCtr,:) = filtfilt(deltafilter.tf.num, 1 , rawDataBySessionNeural.lfpData(ch,:));
                tdbratio.filtdata(chCtr,:) = theta.filtdata(chCtr,:) ./ (delta.filtdata(chCtr,:)+beta.filtdata(chCtr,:));
                tdbratio.filtdata(chCtr,:) = smoothvect(tdbratio.filtdata(chCtr,:), smoothing_kernel_tdb);%smooth data
                tdbratio.baseline(chCtr) = mean(tdbratio.filtdata(chCtr,:));
            end%TDB ratio

            %find all outlier indicies for this channel
            thisChOutInd = [];
            if isfield(rawDataBySessionNeural, 'lfpOutlierInd')%DC added since skipping this step in preprocessing
                for outInd = 1:size(rawDataBySessionNeural.lfpOutlierInd,1)
                    thisChOutInd = [thisChOutInd rawDataBySessionNeural.lfpOutlierInd(outInd,1):rawDataBySessionNeural.lfpOutlierInd(outInd,2)];
                end
            end
            %loop through ripples
            for r = 1:length(ripplesGood(chCtr).midind)

                %threshold for ratio of high frequencies (150-250/250-350 Hz)
                if params.ripple.applyHighFreqRatio
                    % compute psd
                    [Pxx,F] = pwelch(detrend(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))),[],[],[],ripple.samprate);
                    % compute max & mean power in 150-250 and 250-350
                    subfreqs{1} = find(F>=params.ripple.freqNumerator(1) & F<=params.ripple.freqNumerator(2)); %ripple band
                    subfreqs{2} = find(F>params.ripple.freqDenominator(1) & F<=params.ripple.freqDenominator(2)); %above ripple, prob movement
                    for f = 1:2 %for 150-250 or 250-400
                        highFreqVals(1,f) = mean(Pxx(subfreqs{f}));
                    end%f
                    highFreqRatio(r) = highFreqVals(1)./highFreqVals(2);
                    excludeForHighFreq(r) = highFreqRatio(r)<params.ripple.ratioThresh;
                end%high freq ratio

                %average theta/(delta+beta) ratio
                if params.ripple.applyTDBRatio
                    meantdb(r) = mean(tdbratio.filtdata(chCtr,ripperiods(r,1):ripperiods(r,2)));%average tdbratio during time around center of rip
                    excludeForMeanTDB(r) = meantdb(r)>tdbratio.baseline(chCtr);
                end%TDB ratio

                %speed
                if params.ripple.applySpeed
                    meanRipSpeed(r) = mean(rawDataBySessionNeural.speed(find(rawDataBySessionNeural.lfpTime>=ripperiods(r,1) & rawDataBySessionNeural.lfpTime<=ripperiods(r,2))));
                    excludeForSpeed(r) = meanRipSpeed(r)>params.speedTh;
                end%speed

                %noise in raw data
                % remove ripple if any index lies within an outlier period
                if any(ismember(ripperiods(r,1):ripperiods(r,2),thisChOutInd))
                    excludeForRawDataNoise(r) = 1;
                end

                % noise deflection in raw data; find differences in raw data values above a threshold
                if any(abs(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) >= ...
                        std(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))*params.ripple.nstdNoise)
                    excludeForRawDataNoise(r) = 1;
                end%noise deflection
                % % %plot criterion
                % % figure
                % % plot(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))
                % % figure; hold on
                % % plot(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))
                % % plot(ones(1,length(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) * mean(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))))
                % % plot(ones(1,length(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) * std(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))*params.ripple.nstdNoise)
                % % plot(ones(1,length(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) * -std(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))*params.ripple.nstdNoise)

                %spikes present
                % note: looks for spikes from clusters across all channels, not just clusters from the current channel (ch) this loop
                if params.ripple.applyMUA
                    for clu = 1:length(rawDataBySessionNeural.apData)
                        if any(rawDataBySessionNeural.apData(clu).spikeInds/params.samprate >= ripperiods(r,1)/params.lfp_samprate_down &...
                                rawDataBySessionNeural.apData(clu).spikeInds/params.samprate <= ripperiods(r,2)/params.lfp_samprate_down)
                            numSpikes = [];
                            numSpikes = length(find(rawDataBySessionNeural.apData(clu).spikeInds/params.samprate >= ripperiods(r,1)/params.lfp_samprate_down &...
                                rawDataBySessionNeural.apData(clu).spikeInds/params.samprate <= ripperiods(r,2)/params.lfp_samprate_down));
                            %save the spikes for each ripple
                            ripSpikes(r,clu,1:numSpikes) = rawDataBySessionNeural.apData(clu).spikeInds(find(rawDataBySessionNeural.apData(clu).spikeInds/params.samprate >= ripperiods(r,1)/params.lfp_samprate_down &...
                                rawDataBySessionNeural.apData(clu).spikeInds/params.samprate <= ripperiods(r,2)/params.lfp_samprate_down));
                        end%find spikes
                    end%clu
                    if sum(~isnan(ripSpikes(r,:,:)),'all') > 0%spikes
                        excludeForNoSpikesPresent(r) = 0;
                    else%no spikes
                        excludeForNoSpikesPresent(r) = 1;
                    end%if sum(~isnan(ripSpikes(r,:,:)),'all') > 0
                else
                    excludeForNoSpikesPresent(r) = 0;
                end%spikes present

            end%r

            %determine ripples to exclude
            ripplesToExclue = excludeForHighFreq | excludeForMeanTDB | excludeForSpeed | excludeForRawDataNoise | excludeForNoSpikesPresent;
            %save start indices of bad ripples to separate struct
            ripplesBad(chCtr).midind = ripplesGood(chCtr).midind(ripplesToExclue);
            ripplesBad(chCtr).exclusionReason = [excludeForHighFreq(ripplesToExclue); excludeForMeanTDB(ripplesToExclue); ...
                excludeForSpeed(ripplesToExclue);excludeForRawDataNoise(ripplesToExclue); excludeForNoSpikesPresent(ripplesToExclue)];
            %remove ripples that meet exclusion criteria
            ripplesGood(chCtr).startind(ripplesToExclue) = [];
            ripplesGood(chCtr).midind(ripplesToExclue) = [];
            ripplesGood(chCtr).endind(ripplesToExclue) = [];
            ripplesGood(chCtr).starttime(ripplesToExclue) = [];
            ripplesGood(chCtr).midtime(ripplesToExclue) = [];
            ripplesGood(chCtr).endtime(ripplesToExclue) = [];
            ripplesGood(chCtr).peak(ripplesToExclue) = [];
            ripplesGood(chCtr).energy(ripplesToExclue) = [];
            ripplesGood(chCtr).maxthresh(ripplesToExclue) = [];
            %save spikes for this channel
            ripplesGood(chCtr).ripSpikes = ripSpikes(~ripplesToExclue,:,:);

        else

            excludeForHighFreq = [];
            excludeForMeanTDB = [];
            excludeForSpeed = [];
            excludeForRawDataNoise = [];
            excludeForNoSpikesPresent = [];
        end%if any(ripples(ch).midind)

    end%ch

end%apply ripple criteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% add to sturct and save data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawDataBySessionNeural.ripplesGood = ripplesGood;
rawDataBySessionNeural.ripplesBad = ripplesBad;
filename = [saveNeuralPath '\' 'rawDataBySessionNeural.mat'];
save(filename, 'rawDataBySessionNeural', '-v7.3')

%% extract ripples across laps %%
%Note: Rest sessions do not have lap structs

if isfile([saveNeuralPath '\rawDataByLapNeural.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% load session data %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([saveNeuralPath '\rawDataByLapNeural.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% collect ripple info per lap %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:length(rawDataByLapNeural)
        for ch = 1:length(rawDataBySessionNeural.ripplesGood)
            ripLapCtr = 0;
            for r = 1:length(rawDataBySessionNeural.ripplesGood(ch).startind)
                if ismember(rawDataBySessionNeural.ripplesGood(ch).startind(r), rawDataByLapNeural(ii).lfpTime(1):rawDataByLapNeural(ii).lfpTime(end))
                    ripLapCtr = ripLapCtr + 1;
                    rawDataByLapNeural(ii).ripplesGood(ch).startind(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).startind(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).midind(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).midind(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).endind(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).endind(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).starttime(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).starttime(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).midtime(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).midtime(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).endtime(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).endtime(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).peak(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).peak(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).energy(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).energy(r);
                    rawDataByLapNeural(ii).ripplesGood(ch).maxthresh(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).maxthresh(r);
                end%if ismember
            end%r

            if ripLapCtr > 0
                rawDataByLapNeural(ii).ripplesGood(ch).timeRange = rawDataBySessionNeural.ripplesGood(ch).timeRange;
                rawDataByLapNeural(ii).ripplesGood(ch).samprate = rawDataBySessionNeural.ripplesGood(ch).samprate;
                rawDataByLapNeural(ii).ripplesGood(ch).envThreshold = rawDataBySessionNeural.ripplesGood(ch).envThreshold;
                rawDataByLapNeural(ii).ripplesGood(ch).envBaseline = rawDataBySessionNeural.ripplesGood(ch).envBaseline;
                rawDataByLapNeural(ii).ripplesGood(ch).envStd =  rawDataBySessionNeural.ripplesGood(ch).envStd;
                rawDataByLapNeural(ii).ripplesGood(ch).minimum_duration = rawDataBySessionNeural.ripplesGood(ch).minimum_duration;
                rawDataByLapNeural(ii).ripplesGood(ch).chan = rawDataBySessionNeural.ripplesGood(ch).chan;
                rawDataByLapNeural(ii).ripplesGood(ch).ripSpikes = rawDataBySessionNeural.ripplesGood(ch).ripSpikes;
            end%~isempty

        end%ch
    end%lap

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% save data %%%%%
    %%%%%%%%%%%%%%%%%%%%
    filename = [saveNeuralPath '\' 'rawDataByLapNeural.mat'];
    save(filename, 'rawDataByLapNeural', '-v7.3')

end%if isfile([saveNeuralPath '\rawDataByLapNeural.mat'])

%% extract ripples across trials %%
%Note: Rest sessions and active sessions with <= 1 trial do not have trial structs

if isfile([saveNeuralPath '\rawDataByTrialNeural.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% load session data %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([saveNeuralPath '\rawDataByTrialNeural.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% collect ripple info per trial %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for znType = 1:size(rawDataByTrialNeural,1) %1=reward, 2=nonreward, 3 = alt nonreward
        for znNum = 1:size(rawDataByTrialNeural,2)-1
            for tr = 1:size(rawDataByTrialNeural{znType,znNum},2)

                if isempty(rawDataByTrialNeural{znType,znNum})
                    continue
                end

                for ch = 1:length(rawDataBySessionNeural.ripplesGood)
                    ripLapCtr = 0;
                    for r = 1:length(rawDataBySessionNeural.ripplesGood(ch).startind)
                        if ismember(rawDataBySessionNeural.ripplesGood(ch).startind(r), rawDataByTrialNeural{znType,znNum}.lfpTime(1):rawDataByTrialNeural{znType,znNum}.lfpTime(end))
                            ripLapCtr = ripLapCtr + 1;
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).startind(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).startind(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).midind(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).midind(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).endind(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).endind(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).starttime(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).starttime(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).midtime(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).midtime(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).endtime(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).endtime(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).peak(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).peak(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).energy(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).energy(r);
                            rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).maxthresh(ripLapCtr) = rawDataBySessionNeural.ripplesGood(ch).maxthresh(r);
                        end%if ismember
                    end%r

                    if ripLapCtr > 0
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).timeRange = rawDataBySessionNeural.ripplesGood(ch).timeRange;
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).samprate = rawDataBySessionNeural.ripplesGood(ch).samprate;
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).envThreshold = rawDataBySessionNeural.ripplesGood(ch).envThreshold;
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).envBaseline = rawDataBySessionNeural.ripplesGood(ch).envBaseline;
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).envStd =  rawDataBySessionNeural.ripplesGood(ch).envStd;
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).minimum_duration = rawDataBySessionNeural.ripplesGood(ch).minimum_duration;
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).chan = rawDataBySessionNeural.ripplesGood(ch).chan;
                        rawDataByTrialNeural{znType,znNum}.ripplesGood(ch).ripSpikes = rawDataBySessionNeural.ripplesGood(ch).ripSpikes;
                    end%~isempty
                end%ch

            end%tr
        end%znNum
    end%znType

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% save data %%%%%
    %%%%%%%%%%%%%%%%%%%%
    filename = [saveNeuralPath '\' 'rawDataByTrialNeural.mat'];
    save(filename, 'rawDataByTrialNeural', '-v7.3')

end%if isfile([saveNeuralPath '\rawDataByTrialNeural.mat'])

%% plot %%
if plotRipples

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% find channel with max number of ripples %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numGoodRipPerCh = []; numBadRipPerCh = [];
    for ch = 1:length(layerChans)
        numGoodRipPerCh(ch) = length(ripplesGood(ch).midind);
        numBadRipPerCh(ch) = length(ripplesBad(ch).midind);
    end
    [~, maxNumRipChInd] = max(numGoodRipPerCh);

    %plot whole raw data and ripple-filtered traces with stars for start of each detected ripple
    figure; hold on
    chSpace = 500;
    p1 = plot(rawDataBySessionNeural.lfpData(layerChans(maxNumRipChInd),1:length(ripple.filtdata(maxNumRipChInd,:)))*500+1*chSpace);%raw trace
    p2 = plot(ripple.filtdata(maxNumRipChInd,:)*500+1*chSpace);%ripple filtered trace
    p3 = plot(ripple.env(maxNumRipChInd,:)*500+1*chSpace);%ripple env
    p4 = plot(1:length(ripple.filtdata(maxNumRipChInd,:)), ones(1,length(ripple.filtdata(maxNumRipChInd,:)))* ripplesGood(maxNumRipChInd).envBaseline*500+1*chSpace);%env baseline
    p5 = plot(1:length(ripple.filtdata(maxNumRipChInd,:)), ones(1,length(ripple.filtdata(maxNumRipChInd,:)))* ripplesGood(maxNumRipChInd).envThreshold*500+1*chSpace);%env threshold
    for ch = 1:length(layerChans)
        p6 = plot(ripplesGood(ch).midind,ripple.filtdata(maxNumRipChInd,ripplesGood(ch).midind)*500+1*chSpace,'*g');%detected good ripples
        if ~isempty(ripplesBad)
            if ~isempty(ripplesBad(ch).midind)
                p7 = plot(ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(3,:))),...
                    ripple.filtdata(maxNumRipChInd,ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(3,:))))*500+1*chSpace,'*r');%detected speed bad ripples
                p8 = plot(ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(4,:))),...
                    ripple.filtdata(maxNumRipChInd,ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(4,:))))*500+1*chSpace,'*c');%detected raw data bad ripples
                p9 = plot(ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(5,:))),...
                    ripple.filtdata(maxNumRipChInd,ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(5,:))))*500+1*chSpace+10,'*b');%detected MUA bad ripples
            end%if ~isempty(ripplesBad(ch).midind)
        end%if ~isempty(ripplesBad)
    end%ch
    %legend
    legNames = {'RawTrace', 'RipTrace', 'RipEnv', 'EnvBase', 'EnvThresh'};
    if ~isempty(p6)
        legNames{end+1} = 'GoodRip';
    end
    if ~isempty(p7)
        legNames{end+1} = 'SpeedExc';
    end
    if ~isempty(p8)
        legNames{end+1} = 'NoiseExc';
    end
    if ~isempty(p9)
        legNames{end+1} = 'MUAExc';
    end
    leg = legend([p1 p2 p3 p4 p5 p6 p7 p8 p9], legNames,'Box','off',Location='southeast');
    legTitle = get(leg,'title');
    set(legTitle, 'String','Max Channel:', 'FontSize',8);

    %plot spiking on top of filtered ripple trace
    ch = maxNumRipChInd;
    for r = 1:length(ripplesGood(ch).startind)
        %gather spikes for all clusters
        spikeTrainLFPTime = round(squeeze(ripplesGood(ch).ripSpikes(r,:,:)) / params.samprate * params.lfp_samprate_down);

        %initialize
        thisRipPossSp = repmat(1:ripplesGood(ch).endind(r),size(spikeTrainLFPTime,1),1);%matrix for all possible spike times
        thisRipRealSp = zeros(size(thisRipPossSp,1), size(thisRipPossSp,2));%matrix of zeroes for all real spike times

        %find spike times for each cluster
        for clu = 1:size(spikeTrainLFPTime,1)
            thisCluSpInd = ~isnan(spikeTrainLFPTime(clu,:));
            thisRipRealSp(clu,spikeTrainLFPTime(clu,thisCluSpInd)) = ripple.filtdata(ch,spikeTrainLFPTime(clu,thisCluSpInd));
        end
        thisRipRealSp(thisRipRealSp==0) = nan;%turn zeros into NaNs

        %(optional) reduce to only putative pyramidal cells roughly in hippocampus
        hipChansOnly = 1;
        if hipChansOnly
            for clu = 1:size(spikeTrainLFPTime,1)
                if ismember(rawDataBySessionNeural.apData(clu).maxChan,[1:150]) &&...
                        contains(rawDataBySessionNeural.apData(clu).putativeCellType, 'Pyr')
                else
                    thisRipRealSp(clu,:) = nan;
                end%if
            end%chan
        end%pyrLayerChansOnly

        %plot
        colors = cbrewer('div','RdBu', size(spikeTrainLFPTime,1)+3); colors = flip(colors);
        figure; hold on
        plot(thisRipPossSp(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)), rawDataBySessionNeural.lfpData(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)));%raw trace
        plot(thisRipPossSp(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)), ripple.filtdata(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)));%filtered ripple trace
        img = arrayfun( @(x) plot(thisRipPossSp(x,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)),thisRipRealSp(x,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)),'|',...
            'Color', abs(colors(x,:)), 'DisplayName', num2str(x), 'LineWidth', 5), 1:size(spikeTrainLFPTime,1));%each cluster individual color
    end%ripples

end%if plotRipples

end%function
