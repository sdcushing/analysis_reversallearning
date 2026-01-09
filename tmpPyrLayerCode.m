function sessionPyrLayerInfo = getPyrLayerInfo(subj, sessDate, sessNum, params,neuralRawDataPath, saveNeuralPath,  plotPyrLayer)
%last checked JLK 10/15/25

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
    if sum(rawDataBySessionNeural.speedSmooth(randInds(t):randInds(t)+numsecs/0.02-1)) == 0 %samprate ~= 0.02 secs
        tCtr = tCtr + 1;
        randImTime(tCtr) = randInds(t);
    end
end
if length(randImTime)>maxNumEvents; randImTime = randImTime(1:maxNumEvents);end%cut to max number of events


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
freq(4).name = 'spiking';
freq(4).range = [150 240];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find power for all channels to use as min/max values %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find power for each channel
mych = 1:size(rawDataBySessionNeural.lfpData,1);%all channels
all_ch_power_sp = [];
chctr = 0;
all_ch_lfp = [];
all_ch_power = [];

for ch = mych
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
end%ch

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

numChans = 10;
allWSpPowerAmp = []; allWPyrCount = []; allWSpAmpMed = [];
for w = 1:length(all_ch_power_sp)-numChans+1
    thisW = w:w+(numChans-1);
    allWSpPowerAmp(w) = median(all_ch_power_sp(thisW));
    allWPyrCount(w) = length(allPyrCluID(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
    allWSpAmpMed(w) = median(allPyrSpAmpsMn(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
end
allWPyrCount(allWPyrCount==0) = nan;
allWSpAmpMed(allWSpAmpMed==0) = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% loop through windows of channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpPeaks = findpeaks(smooth(allWSpPowerAmp,5),median(smooth(allWSpPowerAmp,5))+std(smooth(allWSpPowerAmp,5)));
tmpPeaks = tmpPeaks.loc;
tmpLay = [];
[~, tmpLay] = max(allWSpPowerAmp(tmpPeaks));
tmpPeaksSansMax = tmpPeaks;
tmpPeaksSansMax(tmpLay) = 1;
[~, tmpLay(2)] = max(allWSpPowerAmp(tmpPeaksSansMax));
pyrLayerCA1 = max(tmpPeaks(tmpLay));
pyrLayerCA3 = min(tmpPeaks(tmpLay));

%%%%%%%%%%%%%%%%%%%%%
%%%%% save data %%%%%
%%%%%%%%%%%%%%%%%%%%%
sessionPyrLayerInfo.pyrLayerCA1 = pyrLayerCA1-10:pyrLayerCA1+9;
sessionPyrLayerInfo.pyrLayerCA3 = pyrLayerCA3-10:pyrLayerCA3+9;
filename = [saveNeuralPath '\' 'sessionPyrLayerInfo.mat'];
save(filename, 'sessionPyrLayerInfo', '-v7.3')

%% Plot %%
if plotPyrLayer

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% plot depth plots of raw + filtered traces %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% plot spike power amplitude %%%%%
    %Note: peaks in smoothed data used to determine pyramidal layers
    figure; hold on; plot(allWSpPowerAmp, '-b'); plot(tmpPeaks(tmpLay), allWSpPowerAmp(tmpPeaks(tmpLay)), '*r');
    plot(smooth(allWSpPowerAmp,5), '-r')
    title(['Pyramidal Layer Identification: ' subj '_' sessDate '_' sessNum], 'Interpreter','none')

    %%%%% find ripple amplitudes by channel %%%%%
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

    %Note: Can add black lines to investigate phase across channels shift
    %plot raw data
    chans = 1:180;
    ytickspacing = [];
    thisch = 0;
    figure; hold on; for ch = chans; thisch = thisch+1; plot(rawDataBySessionNeural.lfpData(ch,randImTime:randImTime+chparams.Fs*numsecs-1)+ch/10'); ytickspacing(thisch) = ch/10;end%can change devisor for spacing
    title('CA1 Raw Traces of Random Immobile Time')
    xticks([1:chparams.Fs-1:chparams.Fs*numsecs-1])
    xticklabels([0 1 2 3])
    xlabel('Time (s)')
    ylabel('Channel')
    yticks([ytickspacing])
    yticklabels([chans])

    %plot filtered data
    plot_spacer = 1000;
    figure; hold on; for ch = 1:length(chans); plot(filtdata(chans(ch),:)'+ch/plot_spacer);end
    title('CA1 Filtered Traces of Random Immobile Time')
    xlabel('Time (s)')
    xticks([1:chparams.Fs-1:chparams.Fs*numsecs-1])
    xticklabels([0 1 2 3])
    ylabel('Channel')
    yticks([1:size(data,2)]/plot_spacer)
    yticklabels([chans])

end%if plotPyrLayer

end%function



