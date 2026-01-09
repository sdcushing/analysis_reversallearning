function sessionPyrLayerInfo = getPyrLayerInfo_temptemp(subj, sessDate, sessNum, params,neuralRawDataPath, saveNeuralPath,  plotPyrLayer)

%% load session data %%
load([neuralRawDataPath '\kilosort4\clusters_allrec.mat'])
load([saveNeuralPath '\rawDataBySessionNeural.mat'])

%% determine random time based on mobility %%
numsecs = 2;
maxNumEvents = 100;
randMobTime = [];
randInds = randi(length(rawDataBySessionNeural.vrTime)-numsecs/0.02-1,100000,1);
tCtr = 0;
for t = 1:length(randInds)
    if min(rawDataBySessionNeural.speedSmooth(randInds(t):randInds(t)+numsecs/0.02-1)) >= 1%0.3
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find power for all channels to use as min/max values %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find power for each channel
mych = 1:size(rawDataBySessionNeural.lfpData,1);%all channels

all_ch_power_th = [];
chctr = 0;
all_ch_lfp = [];
all_ch_power = [];

for ch = mych
    chctr = chctr + 1;
    %get data for this channel
    tmp_lfp = [];
    for tMob = 1:maxNumEvents
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

%find median theta power amplitude, median spike amplitude, and number of putative pyramidal cells for each window of n channels
numChans = 10;
allWThPowerAmp = []; allWPyrCount = []; allWSpAmpMed = []; allWPyrCountPeaks = []; pyrLayerCA1 = []; pyrLayerCA3 = [];
for w = 1:length(all_ch_power_th)-numChans+1
    thisW = w:w+(numChans-1);
    allWThPowerAmp(w) = median(all_ch_power_th(thisW));
    allWPyrCount(w) = length(allPyrCluID(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
    allWSpAmpMed(w) = median(allPyrSpAmpsMn(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
end
allWPyrCount(allWPyrCount==0) = nan;
allWSpAmpMed(allWSpAmpMed==0) = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% identify layers %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%determine windows with ripple amplitude, spike amplitude, and number of pyramidal cells above 2 std
pyrLayers = []; pyrLayerCA1tmp = []; pyrLayerCA3tmp = []; pyrLayerCA1 = []; pyrLayerCA3 = [];
pyrLayers = find(allWPyrCount>=nanmedian(allWPyrCount)+(nanstd(allWPyrCount)*2) & ...
    allWSpAmpMed>=nanmedian(allWSpAmpMed)+(nanstd(allWSpAmpMed)*2));

%plot
figure
hold on
%plot theta
plot(allWThPowerAmp, '-b')%all windows
%plot pyramidal cell counts
plot(allWPyrCount, '-g')%all windows
%plot pyramidal spike amplitude
plot(allWSpAmpMed, '-c')%all windows
%color in CA1 and CA3 peaks
plot(pyrLayers,allWPyrCount(pyrLayers),'or','MarkerFaceColor','r')
plot(pyrLayers,allWSpAmpMed(pyrLayers),'or','MarkerFaceColor','r')
title(['Pyramidal Layer Identification: ' subj '_' sessDate '_' sessNum], 'Interpreter','none')



%% determine random time based on immobility %%
numsecs = 2;
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

for f = 1:size(filtdata,1)
    tmpVal(f) = max(filtdata(f,:)) -  median(filtdata(f,:),2);
end

[~,tmpLay] = max(tmpVal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% plot depth plots of raw + filtered traces %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note: Can add black lines to investigate phase across channels shift
%plot raw data
%chans = pyrLayers-20:2:pyrLayers+20;%by 4 = 40 microns, so -20 to +20 channels = 600 microns
chans = tmpLay-10:tmpLay+10;
%chans = 1:15:160;
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


% % %%%%%%% OLD %%%%%%%%%
% % %find window with peak theta power amplitude
% % [allWThPowerMaxVal, allWThPowerMaxInd] = max(allWThPowerAmp);
% % %find smoothed window with peak pyrmidal cells count and amplitude
% % allWPyrCountSmooth = smooth(allWPyrCount,7);
% % allWPyrCountPeaks = findpeaks(allWPyrCountSmooth,std(allWPyrCount)*2);
% % allWPyrCountPeaks = allWPyrCountPeaks.loc;
% % allWSpAmpMedAbvThresh = find(allWSpAmpMed>nanmedian(allWSpAmpMed)+(nanstd(allWSpAmpMed)*2));
% % %determine channels for CA1 and CA3 layers
% % pyrLayerCA1 = allWPyrCountPeaks(find(allWPyrCountPeaks>allWThPowerMaxInd,1,'first'));
% % possibleCA3Chans = allWSpAmpMedAbvThresh(ismember(allWSpAmpMedAbvThresh,allWPyrCountPeaks));
% % pyrLayerCA3 =  possibleCA3Chans(find(possibleCA3Chans<allWThPowerMaxInd,1,'last'));
% % 
% % %plot
% % figure
% % hold on
% % %plot theta power
% % plot(allWThPowerAmp, '-b')%all windows
% % plot(allWThPowerMaxInd, allWThPowerMaxVal, '*k')%peak window
% % %plot pyramidal cell counts
% % plot(allWPyrCount, '-g')%all windows
% % plot(allWPyrCountSmooth,'Color',[.5 .5 .5])%smoothed windows
% % %plot pyramidal spike amplitude
% % plot(allWSpAmpMed, '-c')%all windows
% % %color in CA1 and CA3 peaks
% % plot(pyrLayerCA1,allWThPowerAmp(pyrLayerCA1),'or','MarkerFaceColor','r')
% % plot(pyrLayerCA1,allWPyrCount(pyrLayerCA1),'or','MarkerFaceColor','r')
% % plot(pyrLayerCA1,allWSpAmpMed(pyrLayerCA1),'or','MarkerFaceColor','r')
% % plot(pyrLayerCA3,allWThPowerAmp(pyrLayerCA3),'or','MarkerFaceColor','r')
% % plot(pyrLayerCA3,allWPyrCount(pyrLayerCA3),'or','MarkerFaceColor','r')
% % plot(pyrLayerCA3,allWSpAmpMed(pyrLayerCA3),'or','MarkerFaceColor','r')
% % title(['Pyramidal Layer Identification: ' subj '_' sessDate '_' sessNum], 'Interpreter','none')

end%function