function [rawDataBySessionNeural] = getRipples(dirs, params, saveNeuralPath, plotRipples)
%adapted from filtereeg2_Intan.m, extractripples3.m,
% ripplepostfileprocess2, findpowerratioripplevsabove2
%last check JLK 10/29/25

%% load session data %%
load([saveNeuralPath '\rawDataBySessionNeural.mat'])
load([saveNeuralPath '\sessionPyrLayerInfo.mat'])

%% Extract Ripples %%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find good ripples per channel %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ch = layerChans(1)
    chCtr = chCtr+1;
    % ripple filter (150-250 Hz)
    ripple.filtdata(chCtr,:) = filtfilt(ripplefilter.kernel, 1 , rawDataBySessionNeural.lfpData(ch,:));
   
    %Hilbert transform ripple data
    hdata = hilbert(ripple.filtdata(chCtr,:));
    ripple.phase(chCtr,:) = angle(hdata);
    ripple.env(chCtr,:) = abs(hdata);
    ripple.env(chCtr,:) = smoothvect(ripple.env(chCtr,:), smoothing_kernel_rip);%smooth env
    clear hdata

    %find baseline, std, and thresh for this channel's ripple env
    envBaseline = mean(ripple.env(chCtr,:));
    envStd = std(ripple.env(chCtr,:));
    envThresh = envBaseline + params.ripple.nstdEnv(end) * envStd;

    %extract the events if this is a valid trace
    if (envThresh > 0) && any(find(ripple.env(chCtr,:)<envBaseline))
        %find ripple events
        tmprip = [];
        tmprip = extractevents(ripple.env(chCtr,:), envThresh, envBaseline, 0, minRipDur, 0)';
        %start, middle (of energy), and end indices and times
        ripples(chCtr).startind = tmprip(:,1);
        ripples(chCtr).midind = tmprip(:,8);
        ripples(chCtr).endind = tmprip(:,2);
        ripples(chCtr).starttime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).startind / ripple.samprate;
        ripples(chCtr).midtime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).midind / ripple.samprate;
        ripples(chCtr).endtime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).endind / ripple.samprate;
        %ripple characteristics
        ripples(chCtr).peak = tmprip(:,3);
        ripples(chCtr).energy = tmprip(:,7);
        ripples(chCtr).maxthresh = (tmprip(:,9) - envBaseline) / envStd;
    else
        ripples(chCtr).startind = [];
        ripples(chCtr).midind = [];
        ripples(chCtr).endind = [];
        ripples(chCtr).starttime = [];
        ripples(chCtr).midtime = [];
        ripples(chCtr).endtime = [];
        ripples(chCtr).peak = [];
        ripples(chCtr).energy = [];
        ripples(chCtr).maxthresh = [];
    end%> thresh + baseline

    %add other info to ripples struct
    ripples(chCtr).timeRange = [0 length(ripple.env(chCtr,:))/ripple.samprate] + rawDataBySessionNeural.lfpTime(1);
    ripples(chCtr).samprate = ripple.samprate;
    ripples(chCtr).envThreshold = envThresh;
    ripples(chCtr).envBaseline = envBaseline;
    ripples(chCtr).envStd = envStd;
    ripples(chCtr).minimum_duration = params.ripple.minRipDur;
    ripples(chCtr).chan = ch;

    %find ripples that meet criteria
    if any(ripples(chCtr).midind) && params.ripple.applyCriteria
        %find time around ripple middle (2 s total)
        ripperiods = [ripples(chCtr).midind-params.ripple.timeAroundRip*ripple.samprate, ripples(chCtr).midind+params.ripple.timeAroundRip*ripple.samprate];
        highFreqVals = [];
        highFreqRatio = [];
        meantdb = [];
        excludeForHighFreq = zeros(1,length(ripples(chCtr).midind));
        excludeForMeanTDB = zeros(1,length(ripples(chCtr).midind));
        excludeForSpeed = zeros(1,length(ripples(chCtr).midind));
        excludeForRawDataNoise = zeros(1,length(ripples(chCtr).midind));
        excludeForNoSpikesPresent = zeros(1,length(ripples(chCtr).midind));
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

        %loop through ripples
        for r = 1:length(ripples(chCtr).midind)
           
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

            %noise deflection in raw data; find differences in raw data values above a threshold
            if any(abs(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) >= ...
                    std(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))*params.ripple.nstdNoise)
                excludeForRawDataNoise(r) = 1;
            else
                excludeForRawDataNoise(r) = 0;
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
                        ripSpikes(r,clu,1:numSpikes) = find(rawDataBySessionNeural.apData(clu).spikeInds/params.samprate >= ripperiods(r,1)/params.lfp_samprate_down &...
                            rawDataBySessionNeural.apData(clu).spikeInds/params.samprate <= ripperiods(r,2)/params.lfp_samprate_down);
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
        ripplesBad(chCtr).startind = ripples(chCtr).startind(ripplesToExclue);
        ripplesBad(chCtr).exclusionReason = [excludeForHighFreq(ripplesToExclue); excludeForMeanTDB(ripplesToExclue); ...
            excludeForSpeed(ripplesToExclue);excludeForRawDataNoise(ripplesToExclue); excludeForNoSpikesPresent(ripplesToExclue)];
        %remove ripples that meet exclusion criteria
        ripples(chCtr).startind(ripplesToExclue) = [];
        ripples(chCtr).midind(ripplesToExclue) = [];
        ripples(chCtr).endind(ripplesToExclue) = [];
        ripples(chCtr).starttime(ripplesToExclue) = [];
        ripples(chCtr).midtime(ripplesToExclue) = [];
        ripples(chCtr).endtime(ripplesToExclue) = [];
        ripples(chCtr).peak(ripplesToExclue) = [];
        ripples(chCtr).energy(ripplesToExclue) = [];
        ripples(chCtr).maxthresh(ripplesToExclue) = [];
        %save spikes for this channel
        ripples(chCtr).ripSpikes = ripSpikes(~ripplesToExclue,:,:);

    else
        excludeForHighFreq = [];
        excludeForMeanTDB = [];
        excludeForSpeed = [];
        excludeForRawDataNoise = [];
        excludeForNoSpikesPresent = [];
    end%if any(ripples(ch).midind)

end%ch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find channel with max number of ripples %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numRipPerCh = [];
for ch = 1:length(layerChans)
    numRipPerCh(ch) = length(ripples(ch).startind);
end
[maxNumRipChVal, maxNumRipChInd] = max(numRipPerCh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find all unique ripples across channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

allRipStInd = []; allRipCh = [];
for ch = 1:length(layerChans)
    allRipStInd = [allRipStInd; ripples(ch).startind];
end%ch
[allRipStInd, I] = sort(allRipStInd);
if length(allRipStInd) ~= sum(numRipPerCh)
    sprintf('Wrong number of ripples!')
end%if

allUniqRipStInd = allRipStInd(diff(allRipStInd)>(params.ripple.minRipDur * ripples(maxNumRipChInd).samprate * 4));


%% plot %%
if plotRipples
    tmpInd = [];
    tmpRipInd = NaN(length(layerChans),length(maxNumRipChVal));
    tmpRipCh = NaN(length(layerChans),length(maxNumRipChVal));
    
    %find all channels with a ripple around the same time
    for ch = 1:length(layerChans)
        for i = 1:maxNumRipChVal
            tmpInd = find(ripples(ch).startind >= (ripples(maxNumRipChInd).startind(i) - (params.ripple.minRipDur * ripples(maxNumRipChInd).samprate)) & ...
                ripples(ch).startind <= (ripples(maxNumRipChInd).startind(i) + (params.ripple.minRipDur * ripples(maxNumRipChInd).samprate)),1,'first');
            if ~isempty(tmpInd)
                tmpRipInd(ch,i) = tmpInd;
                tmpRipCh(ch,i) = ch;
            end
        end%i
    end%ch
    
    %plot ripple found across multiple channels
    figure; hold on
    i = maxNumRipChVal;%ripple index to plot
    for chan = 1:length(layerChans)
        if ~isempty(ripples(chan).startind) && tmpRipInd(chan,i)~=0
            plot(rawDataBySessionNeural.lfpData(layerChans(chan),ripples(maxNumRipChInd).startind(i)-1000:ripples(maxNumRipChInd).startind(i)+1000)*500+chan*50)
            plot(ripple.filtdata(chan,ripples(maxNumRipChInd).startind(i)-1000:ripples(maxNumRipChInd).startind(i)+1000)*500+chan*50)
        end
    end
    
    %plot all channels regardless of ripple detected
    figure; hold on
    for chan = 1:length(layerChans)
            plot(rawDataBySessionNeural.lfpData(layerChans(chan),ripples(maxNumRipChInd).startind(i)-1000:ripples(maxNumRipChInd).startind(i)+1000)*500+chan*50)
            plot(ripple.filtdata(chan,ripples(maxNumRipChInd).startind(i)-1000:ripples(maxNumRipChInd).startind(i)+1000)*500+chan*50)
    end

    %plot whole raw data and ripple-filtered traces with stars for start of each detected ripple
    figure; hold on
    %max chan
    chSpace = 500;
    p1 = plot(rawDataBySessionNeural.lfpData(layerChans(maxNumRipChInd),:)*500+1*chSpace);%raw trace
    p2 = plot(ripple.filtdata(maxNumRipChInd,:)*500+1*chSpace);%ripple filtered trace
    p3 = plot(ripple.env(maxNumRipChInd,:)*500+1*chSpace);%ripple env
    p4 = plot(allUniqRipStInd,ripple.filtdata(maxNumRipChInd,allUniqRipStInd)*500+1*chSpace,'*y');%detected good ripples across all channels
    p5 = plot(ripples(maxNumRipChInd).startind,ripple.filtdata(maxNumRipChInd,ripples(maxNumRipChInd).startind)*500+1*chSpace,'*g');%detected good ripples
    p6 = plot(ripplesBad(maxNumRipChInd).startind(logical(ripplesBad(maxNumRipChInd).exclusionReason(3,:))),...
        ripple.filtdata(maxNumRipChInd,ripplesBad(maxNumRipChInd).startind(logical(ripplesBad(maxNumRipChInd).exclusionReason(3,:))))*500+1*chSpace,'*r');%detected speed bad ripples
    p7 = plot(ripplesBad(maxNumRipChInd).startind(logical(ripplesBad(maxNumRipChInd).exclusionReason(4,:))),...
        ripple.filtdata(maxNumRipChInd,ripplesBad(maxNumRipChInd).startind(logical(ripplesBad(maxNumRipChInd).exclusionReason(4,:))))*500+1*chSpace,'*c');%detected raw data bad ripples
    p8 = plot(ripplesBad(maxNumRipChInd).startind(logical(ripplesBad(maxNumRipChInd).exclusionReason(5,:))),...
        ripple.filtdata(maxNumRipChInd,ripplesBad(maxNumRipChInd).startind(logical(ripplesBad(maxNumRipChInd).exclusionReason(5,:))))*500+1*chSpace+10,'*b');%detected MUA bad ripples
    p9 = plot(1:length(ripple.filtdata(maxNumRipChInd,:)), ones(1,length(ripple.filtdata(maxNumRipChInd,:)))* ripples(maxNumRipChInd).envBaseline*500+1*chSpace);%env baseline
    p10 = plot(1:length(ripple.filtdata(maxNumRipChInd,:)), ones(1,length(ripple.filtdata(maxNumRipChInd,:)))* ripples(maxNumRipChInd).envThreshold*500+1*chSpace);%env threshold
    %alt chan
    altRipChInd = 19;
    p11 = plot(rawDataBySessionNeural.lfpData(layerChans(altRipChInd),:)*500+2*chSpace);%raw trace
    p12 = plot(ripple.filtdata(altRipChInd,:)*500+2*chSpace);%ripple filtered trace
    p13 = plot(ripple.env(altRipChInd,:)*500+2*chSpace);%ripple env
    p14 = plot(ripples(altRipChInd).startind,ripple.filtdata(altRipChInd,ripples(altRipChInd).startind)*500+2*chSpace,'*g');%detected good ripples
    p15 = plot(ripplesBad(altRipChInd).startind(logical(ripplesBad(altRipChInd).exclusionReason(2,:))),...
        ripple.filtdata(altRipChInd,ripplesBad(altRipChInd).startind(logical(ripplesBad(altRipChInd).exclusionReason(2,:))))*500+2*chSpace,'*r');%detected TDB bad ripples
    p16 = plot(ripplesBad(altRipChInd).startind(logical(ripplesBad(altRipChInd).exclusionReason(3,:))),...
        ripple.filtdata(altRipChInd,ripplesBad(altRipChInd).startind(logical(ripplesBad(altRipChInd).exclusionReason(3,:))))*500+2*chSpace,'*c');%detected raw data bad ripples
    p17 = plot(ripplesBad(altRipChInd).startind(logical(ripplesBad(altRipChInd).exclusionReason(4,:))),...
        ripple.filtdata(altRipChInd,ripplesBad(altRipChInd).startind(logical(ripplesBad(altRipChInd).exclusionReason(4,:))))*500+2*chSpace+10,'*b');%detected MUA bad ripples
    p18 = plot(1:length(ripple.filtdata(altRipChInd,:)), ones(1,length(ripple.filtdata(altRipChInd,:)))* ripples(altRipChInd).envBaseline*500+2*chSpace);%env baseline
    p19 = plot(1:length(ripple.filtdata(altRipChInd,:)), ones(1,length(ripple.filtdata(altRipChInd,:)))* ripples(altRipChInd).envThreshold*500+2*chSpace);%env threshold
    %legend
    legNames = {'RawTrace', 'RipTrace', 'RipEnv', 'GoodRipAllCh', 'GoodRipThisCh'};
    if ~isempty(p6)
        legNames{end+1} = 'SpeedExc';
    end
    if ~isempty(p7)
        legNames{end+1} = 'NoiseExc';
    end
    if ~isempty(p8)
        legNames{end+1} = 'MUAExc';
    end
    legNames{end+1} = 'EnvBase'; legNames{end+1} = 'EnvThresh';
    leg = legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10], legNames,'Box','off',Location='southeast');
    legTitle = get(leg,'title');
    set(legTitle, 'String','Max Channel:', 'FontSize',8);

end%if plotRipples

%% save %%
if maxNumRipChVal>0
    rawDataBySessionNeural.ripples = ripples(maxNumRipChInd);
else
    rawDataBySessionNeural.ripples = [];
end
filename = [saveNeuralPath '\' 'rawDataBySessionNeural.mat'];
save(filename, 'rawDataBySessionNeural', '-v7.3')

end%function
