function sessionPyrLayerInfo = getPyrLayerInfo(subj, sessDate, sessNum, params,neuralRawDataPath, saveNeuralPath,  plotPyrLayer, selectManually)
%last checked 10/15/25 JLK

%% load session data %%
load([neuralRawDataPath '\kilosort4\clusters_allrec.mat'])
load([saveNeuralPath '\rawDataBySessionNeural.mat'])

%% determine random time based on immobility %%
numsecs = 1;
maxNumEvents = 100;
randImTime = [];
randInds = randi(length(rawDataBySessionNeural.vrTime)-numsecs/0.02-1,100000,1);
tCtr = 0;
for t = 1:length(randInds)
    if max(rawDataBySessionNeural.speedSmooth(randInds(t):randInds(t)+numsecs/0.02-1)) < 0.3 %samprate ~= 0.02 secs; moving = above 0.3
        %if sum(rawDataBySessionNeural.speedSmooth(randInds(t):randInds(t)+numsecs/0.02-1)) == 0 %samprate ~= 0.02 secs
        tCtr = tCtr + 1;
        randImTime(tCtr) = randInds(t);
    end
end
if length(randImTime)>maxNumEvents; randImTime = randImTime(1:maxNumEvents);end%cut to max number of events

%% determine random time based on mobility %%
numsecs = 2;
maxNumEvents = 100;
randMobTime = [];
randInds = randi(length(rawDataBySessionNeural.vrTime)-numsecs/0.02-1,100000,1);
tCtr = 0;
for t = 1:length(randInds)
    if min(rawDataBySessionNeural.speedSmooth(randInds(t):randInds(t)+numsecs/0.02-1)) >= 0.3 %samprate ~= 0.02 secs; moving = above 0.3
        tCtr = tCtr + 1;
        randMobTime(tCtr) = randInds(t);
    end
end
if length(randMobTime)>maxNumEvents; randMobTime = randMobTime(1:maxNumEvents);end%cut to max number of events

%% calculate power for each channel %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% chronux parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set chronux params
chparams.fpass = [1 400];
chparams.tapers = [2 3];
chparams.Fs = params.lfp_samprate;
chparams.pad = 0;
chparams.trialave = 0;

%determine frequency ranges to look at
%freq 1, 2, 3 would be theta, slow gamma, fast gamma, respectively
freq = [];
freq(1).name = 'theta';
freq(1).range = [6 10];
freq(1).xlim = [1 2];

freq(4).name = 'spiking';
freq(4).range = [150 240];
freq(4).xlim = [3 4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find power for all channels to use as min/max values %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find power for each channel
mych = 1:size(rawDataBySessionNeural.lfpData,1);%all channels
%mych = 1:round(size(rawDataBySessionNeural.lfpData,1)/2);%bottom half of probe (CA3)
%mych = round(size(rawDataBySessionNeural.lfpData,1)/2):str2double(lfp_meta.nSavedChans)-1;%top half of probe (CA1)
all_ch_power_th = [];
all_ch_power_sp = [];
for fr = [1 4]
    chctr = 0;
    all_ch_lfp = [];
    all_ch_power = [];

    for ch = mych
        if fr == 1
            chctr = chctr + 1;
            %get data for this channel
            tmp_lfp = [];
            for tMob = 1:length(randMobTime)
                tmp_lfp(tMob,:) = rawDataBySessionNeural.lfpData(ch,randMobTime(tMob):randMobTime(tMob)+chparams.Fs*numsecs-1);
            end

            %calculate spectra using Chronux
            %S = frequency bins x windows
            [S, f] = mtspectrumc(tmp_lfp',chparams);

            %find median power across windows; tmp_power = frequency bins x 1
            df = 2*(size(tmp_lfp,1))*chparams.tapers(2);%degrees of freedom = 2 * #trials * #tapers
            tmp_power = median(S,2);%use median across events to avoid impact of outlier events

            %save data to all channels struct
            all_ch_lfp(chctr,:,:) = tmp_lfp;
            all_ch_power(chctr,:) = 10*(log10(tmp_power)-psi(df/2)+log(df/2)); %log transformed and bias corrected and X10 bel-->decibel

            whichf = f>=freq(1).range(1) & f<=freq(1).range(2);
            all_ch_power_th(chctr) = mean(all_ch_power(chctr,whichf));
        elseif fr == 4
            chctr = chctr + 1;
            %get data for this channel
            tmp_lfp = [];
            for tImmob = 1:length(randImTime)
                tmp_lfp(tImmob,:) = rawDataBySessionNeural.lfpData(ch,randImTime(tImmob):randImTime(tImmob)+chparams.Fs*numsecs-1);
            end

            %calculate spectra using Chronux
            %S = frequency bins x windows
            [S, f] = mtspectrumc(tmp_lfp',chparams);

            %find median power across windows; tmp_power = frequency bins x 1
            df = 2*(size(tmp_lfp,1))*chparams.tapers(2);%degrees of freedom = 2 * #trials * #tapers
            tmp_power = median(S,2);%use median across events to avoid impact of outlier events

            %save data to all channels struct
            all_ch_lfp(chctr,:,:) = tmp_lfp;
            all_ch_power(chctr,:) = 10*(log10(tmp_power)-psi(df/2)+log(df/2)); %log transformed and bias corrected and X10 bel-->decibel

            whichf = f>=freq(4).range(1) & f<=freq(4).range(2);
            all_ch_power_sp(chctr) = mean(all_ch_power(chctr,whichf));
        end%if fr
    end%ch
end%fr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% collect cluster info %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allPyrSpAmpsMn = [];
allPyrMaxChan = [];
allPyrCluID = [];
allPyrNumSp = [];
allIntCluID = [];
allIntMaxChan = [];
allWideCluID = [];
for clu = 1:length(clusters_allrec)
    if sum([rawDataBySessionNeural.apData.ID]==clusters_allrec(clu).ID)>0
        if contains(rawDataBySessionNeural.apData([rawDataBySessionNeural.apData.ID]==clusters_allrec(clu).ID).putativeCellType, 'Pyr')
            allPyrCluID = [allPyrCluID clusters_allrec(clu).ID];
            allPyrMaxChan = [allPyrMaxChan clusters_allrec(clu).maxChan];
            allPyrSpAmpsMn = [allPyrSpAmpsMn mean(clusters_allrec(clu).spikeAmps)];
            allPyrNumSp = [allPyrNumSp length(clusters_allrec(clu).spikeAmps)];
        elseif contains(rawDataBySessionNeural.apData([rawDataBySessionNeural.apData.ID]==clusters_allrec(clu).ID).putativeCellType, 'Narrow')
            allIntCluID = [allIntCluID clusters_allrec(clu).ID];
            allIntMaxChan = [allIntMaxChan clusters_allrec(clu).maxChan];
        elseif contains(rawDataBySessionNeural.apData([rawDataBySessionNeural.apData.ID]==clusters_allrec(clu).ID).putativeCellType, 'Wide')
            allWideCluID = [allWideCluID clusters_allrec(clu).ID];
        end
    end
end
if (length(allPyrCluID) + length(allIntCluID) + length(allWideCluID)) ~= length(clusters_allrec)
    sprintf('Too many clusters!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% loop through windows of channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find median ripple amplitude, median spike amplitude, and number of putative pyramidal cells for each window of n channels
numChans = 10;
allWThPowerAmp = []; allWSpPowerAmp = []; allWPyrCount = []; allWSpAmpMed = [];
for w = 1:length(all_ch_power_sp)-numChans+1
    thisW = w:w+(numChans-1);
    allWThPowerAmp(w) = median(all_ch_power_th(thisW));
    allWSpPowerAmp(w) = median(all_ch_power_sp(thisW));
    allWPyrCount(w) = length(allPyrCluID(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
    allWSpAmpMed(w) = median(allPyrSpAmpsMn(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
end
allWPyrCount(allWPyrCount==0) = nan;
allWSpAmpMed(allWSpAmpMed==0) = nan;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find layers %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%%%% find peaks in the spike power amplitude %%%%%
smoothVal = 5;
tmpPeaks = findpeaks(smooth(allWSpPowerAmp,smoothVal),median(smooth(allWSpPowerAmp,smoothVal))+std(smooth(allWSpPowerAmp,smoothVal))*0.25);
tmpPeaks = tmpPeaks.loc;

%%%%% find CA3 layer %%%%%
tmpLay = [];
[~, maxThPowerW] = max(allWThPowerAmp);
tmpCA3Peaks = tmpPeaks(abs(tmpPeaks-maxThPowerW)<=15);
[~, tmpLay] = max(allWSpPowerAmp(tmpCA3Peaks));
pyrLayerCA3 = tmpCA3Peaks(tmpLay);

%%%%% find CA1 layer %%%%%
tmpLay = [];
tmpCA1Peaks = tmpPeaks(tmpPeaks>pyrLayerCA3+20);
[~, tmpLay] = max(allWSpPowerAmp(tmpCA1Peaks));
pyrLayerCA1 = tmpCA1Peaks(tmpLay);

%%%%% save data to struct %%%%%
%note: below will override these values if they are empty or selectManually = 1
sessionPyrLayerInfo.pyrLayerCA1 = pyrLayerCA1-10:pyrLayerCA1+9;
sessionPyrLayerInfo.pyrLayerCA3 = pyrLayerCA3-10:pyrLayerCA3+9;

%% Plot %%
if plotPyrLayer || selectManually || isempty(pyrLayerCA3) || isempty(pyrLayerCA1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% plot depth plots of layer detection and LFP traces %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure

    %%%%% plot raw data %%%%%
    %find ripple amplitudes by channel
    filt_freq = [150 240];
    tmpData = [];
    for t = 1:length(randImTime)
        %filt params
        tmpData(t,:,:) = rawDataBySessionNeural.lfpData(:,randImTime(t):randImTime(t)+chparams.Fs*numsecs-1);
    end
    tmpDataMn = squeeze(mean(tmpData,1))';

    %build filt
    N = 3;
    Wn = filt_freq ./ (chparams.Fs/2);
    [b,a] = butter(N,Wn);
    %filter data in spike frequency range
    filtdata = filtfilt(b, a, tmpDataMn);
    filtdata = -1*filtdata; %flip signal so spikes go up
    if size(filtdata,2)>1
        filtdata=filtdata';
    end

    %plot raw data
    subplot(1,2,1); hold on
    chans = 1:140;
    ytickspacing = [];
    thisch = 0;
    for ch = chans; thisch = thisch+1; plot(rawDataBySessionNeural.lfpData(ch,randImTime:randImTime+chparams.Fs*numsecs-1)+ch/10'); ytickspacing(thisch) = ch/10;end%can change devisor for spacing
    title('CA1 Raw Traces of Random Immobile Time')
    xticks([1:chparams.Fs-1:chparams.Fs*numsecs-1])
    xticklabels([0 1 2 3])
    xlabel('Time (s)')
    ylabel('Channel')
    yticks([ytickspacing])
    yticklabels([chans])

    % % uncomment to include filtered lfps across channels
    % % %plot filtered data
    % % plot_spacer = 1000;
    % % figure; hold on; for ch = 1:length(chans); plot(filtdata(chans(ch),:)'+ch/plot_spacer);end
    % % title('CA1 Filtered Traces of Random Immobile Time')
    % % xlabel('Time (s)')
    % % xticks([1:chparams.Fs-1:chparams.Fs*numsecs-1])
    % % xticklabels([0 1 2 3])
    % % ylabel('Channel')
    % % yticks([1:size(filtdata,2)]/plot_spacer)
    % % yticklabels([chan


    %%%%% plot spike power amplitude %%%%%
    %Note: peaks in smoothed data used to determine pyramidal layers
    %plot spike power
    subplot(1,2,2); hold on
    plot(allWSpPowerAmp, '-b'); plot([pyrLayerCA1 pyrLayerCA3], allWSpPowerAmp([pyrLayerCA1 pyrLayerCA3]), '*r'); plot(smooth(allWSpPowerAmp,smoothVal), '-r')
    %plot theta power
    plot(allWThPowerAmp, '-g'); plot(maxThPowerW, allWThPowerAmp(maxThPowerW), '*g')
    title(['Pyramidal Layer Identification: ' subj '_' sessDate '_' sessNum], 'Interpreter','none')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% input pyramidal layer info if empty or manual selection %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(pyrLayerCA3) || isempty(pyrLayerCA1) || selectManually

        % instruct user to click on the unit-count plot (ax2) only
        hFig = gcf;
        try
            targetAx = ax2; % unit count / metric plot
        catch
            % fallback: use current axes
            targetAx = gca;
        end
        figure(hFig); % bring to front
        % make ax2 the current axes so ginput uses its coordinate system
        try
            set(hFig, 'CurrentAxes', targetAx);
        catch
            axes(targetAx);
        end
        title(targetAx, 'Manual selection (ax2): click 4 points in order: CA1 start, CA1 end, CA3 start, CA3 end');
        
        %draw on the plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% XIAO HELP ME!!!! %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        drawnow;

        % collect and draw 4 x-positions with mouse clicks on ax2 only
        [xv,~] = ginput(4);
        xv = round(xv);
        CA1_start_ch = xv(1);
        CA1_end_ch = xv(2);
        CA3_start_ch = xv(3);
        CA3_end_ch = xv(4);
        hold(targetAx, 'on');
        xline(targetAx, CA1_start_ch, '--r', {'CA1 start (manual)'});
        xline(targetAx, CA1_end_ch, '--r', {'CA1 end (manual)'});
        xline(targetAx, CA3_start_ch, '--g', {'CA3 start (manual)'});
        xline(targetAx, CA3_end_ch, '--g', {'CA3 end (manual)'});
        hold(targetAx, 'off');
    end%if isempty(pyrLayerCA3) || isempty(pyrLayerCA1) || selectManually

end%if plotPyrLayer

%%%%%%%%%%%%%%%%%%%%%
%%%%% save data %%%%%
%%%%%%%%%%%%%%%%%%%%%
filename = [saveNeuralPath '\' 'sessionPyrLayerInfo.mat'];
save(filename, 'sessionPyrLayerInfo', '-v7.3')

end%function



