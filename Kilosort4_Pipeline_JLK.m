function Kilosort4_Pipeline_JLK(subj, sessDate, neuralRawDataPath)
%Note: getKilosort4Out.m creates a single AP .bin file for all recordings on a
%given day and runs kilosort4 to produce the kilosort output files used in
%this function. The separate function was created to allow for manually
%spike sorting via phy2 before starting this function.

%% Kilosort4 Pipeline for Neuropixels %%
%NOTE THAT THE KILOSORT SPIKE SORTING PIPELINE COMBINES ALL RECORDINGS
% TOGETHER TO DETECT CLUSTERS THROUGHOUT THE RECORDING DAY, THEN CREATES STRUCTS THAT
% SPLIT THE DATA INTO INDIVIDUAL RECORDING SESSIONS AND THAT COMBINED ALL
% DATA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% collect property information for each recording %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpDirs = dir([neuralRawDataPath]);% '\' sprintf('*%s_%s*',subj,sessDate)
props = [];
for r = 1:size(tmpDirs,1)
    %read ap meta data
    ap_binName = sprintf([tmpDirs(r).name '_t0.imec0.ap.bin']);
    ap_path = [tmpDirs(r).folder '\' tmpDirs(r).name '\' sprintf([tmpDirs(r).name '_imec0'])];

    ap_meta = SGLX_readMeta.ReadMeta(ap_binName, ap_path);

    %read ap data
    samp0 = 0;
    ap_samprate = str2double(ap_meta.imSampRate);%Hz
    nSamp = (str2double(ap_meta.fileTimeSecs)*ap_samprate);

    if strcmp(subj,'JK15') && strcmp(sessDate,'250530') && r == 2%missing data for JK15_250530_3 (session 2)
        nSamp = nSamp-1539684;
    end

    %save info
    props.fileNames = tmpDirs;
    props.recLength(r) = nSamp;
    props.sampRate = ap_samprate;
    props.numChan = str2double(ap_meta.nSavedChans);
    props.fileNums(r) = str2double(tmpDirs(r).name(end-3));
end%r

kilosort_path = [neuralRawDataPath '\' 'kilosort4'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% create raw cluster structures using phy2 output %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adapted from makeClusterStructure.m

%%%%% read clustered information %%%%%
addpath('\\ad.gatech.edu\bme\labs\singer\Danielle\code\CellExplorer\toolboxes\npy-matlb');
spikeInds = readNPY([kilosort_path '\' 'spike_times.npy']);
spikeAmps = readNPY([kilosort_path '\' 'amplitudes.npy']);
spikeID = readNPY([kilosort_path '\' 'spike_clusters.npy']);
[clusterID, clusterGroup] = readClusterGroupsCSV([kilosort_path '\' 'cluster_group.tsv']);
templates = readNPY([kilosort_path '\' 'templates.npy']);
spikeTemplates = readNPY([kilosort_path '\' 'spike_templates.npy']);
channelMap = readNPY([kilosort_path '\' 'channel_map.npy']);

%%%%% only units classified as "good" %%%%%
goodUnits = clusterID(clusterGroup == 2);

%%%%% get templates for each cluster - noise, mua, and good %%%%%
tempPerUnit = findTempForEachClu(spikeID, spikeTemplates);

%check length of tempPerUnit and spikeID are equal
if ~isequal(length(tempPerUnit)-1, clusterID(end))
    error('Length of templates and identified clusters do not match')
end

%%%%% get max channel per cluster based on max template amplitude %%%%%
[~,max_site] = max(max(abs(templates),[],2),[],3);
templateMaxChan = channelMap(max_site); %0-based, template 0 is at ind 1 - max channel of each template
unitMaxChan = templateMaxChan(tempPerUnit(~isnan(tempPerUnit))+1); %only valid templates, +1 because template is 0-based
%unitMaxChan = templateMaxChan(tempPerUnit+1); % +1 because template is 0-based
unitMaxChan = double(unitMaxChan(clusterGroup == 2)); %only good units

%%%%% loop through recordings %%%%%
rawclusters = [];
numShanks = 1;
brainReg = {'CA1+CA3'};
for r = 1:size(props.fileNames,1)
    %create structure
    tmprawclusters = struct('ID', num2cell(goodUnits), ...
        'spikeInds', repmat({[]}, 1, length(goodUnits)),...
        'spikeAmps', repmat({[]}, 1, length(goodUnits)),...
        'sampRate', num2cell(props.sampRate*ones(1, length(goodUnits))), ...
        'maxChan', num2cell(unitMaxChan'), 'info', repmat({'pre quality control metrics'}, 1, length(goodUnits)),...
        'numShanks',num2cell(numShanks*ones(1, length(goodUnits))), 'brainReg', repmat({brainReg}, 1, length(goodUnits)), ...
        'numChan', num2cell(props.numChan*ones(1,length(goodUnits))));

    [tmprawclusters(1:length(goodUnits)).index] = deal(subj);
    [tmprawclusters(1:length(goodUnits)).files] = deal(props.fileNums);

    for clu = 1:length(goodUnits)
        if r == 1
            tempSpikeInds{clu} = spikeInds(spikeID == goodUnits(clu));
            tempSpikeInds{clu} = double(tempSpikeInds{clu});
            tempSpikeAmps{clu} = spikeAmps(spikeID == goodUnits(clu));
            tempSpikeAmps{clu} = double(tempSpikeAmps{clu});
        else
            tempSpikeInds{clu} = tempSpikeInds{clu} - props.recLength(r-1); %align to start time of this recording
            tmprawclusters(clu).spikeInds = [];
        end

        tmprawclusters(clu).spikeInds = tempSpikeInds{clu}(tempSpikeInds{clu} <= props.recLength(r));
        tmprawclusters(clu).spikeAmps = tempSpikeAmps{clu}(tempSpikeInds{clu} <= props.recLength(r));

        if r <  size(props.fileNames,1)
            tempSpikeInds{clu} = tempSpikeInds{clu}(tempSpikeInds{clu} > props.recLength(r));
            tempSpikeAmps{clu} = tempSpikeAmps{clu}(tempSpikeInds{clu} > props.recLength(r));
        end
        tmprawclusters(clu).file = props.fileNums(r);
    end
    %output = rawclusters{session}(cluster)
    rawclusters{r} = tmprawclusters;
end

%make structure with all spike times from all recordings
rawclusters_allrec = struct('ID', num2cell(goodUnits), ...
    'spikeInds', repmat({[]}, 1, length(goodUnits)),...
    'spikeAmps', repmat({[]}, 1, length(goodUnits)),...
    'sampRate', num2cell(props.sampRate*ones(1, length(goodUnits))), ...
    'maxChan', num2cell(unitMaxChan'), 'info', repmat({'all files. pre quality control metrics'}, 1, length(goodUnits)), ...
    'numShanks',num2cell(numShanks*ones(1, length(goodUnits))), 'brainReg', repmat({brainReg}, 1, length(goodUnits)), ...
    'numChan', num2cell(props.numChan*ones(1,length(goodUnits))));
[rawclusters_allrec(1:length(goodUnits)).index] = deal([subj '_' sessDate]);
[rawclusters_allrec(1:length(goodUnits)).files] = deal(props.fileNums);

for clu = 1:length(goodUnits)
    tempSpikeInds{clu} = spikeInds(spikeID == goodUnits(clu));
    tempSpikeInds{clu} = double(tempSpikeInds{clu});
    rawclusters_allrec(clu).spikeInds = tempSpikeInds{clu};
    tempSpikeAmps{clu} = spikeAmps(spikeID == goodUnits(clu));
    tempSpikeAmps{clu} = double(tempSpikeAmps{clu});
    rawclusters_allrec(clu).spikeAmps = tempSpikeAmps{clu};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get waveform information %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adapted from getWaveForms_K2 and subfunctions makeWFstructure_fromRaw.m
% and getstableclustertimes_gauss_K2

%%%%% define time around spike to collect %%%%%
tAroundSpike = [0.001*props.sampRate .002*props.sampRate]; %1ms before and 2ms after
wfWin = -tAroundSpike(1):tAroundSpike(2);
nWFsamps = numel(wfWin);

%%%%% build filter %%%%%
hpFilt = designfilt('highpassiir', 'StopbandFrequency', 100, ...
    'PassbandFrequency', 500, 'StopbandAttenuation', 60, ...
    'PassbandRipple', 1, 'SampleRate', props.sampRate, 'DesignMethod', 'butter');

%%%%% %define some parameters %%%%%
nWFsToLoad = 1000; %number of waveforms to load from the recording
nSamp = sum(props.recLength); %number of total samples for this recording day
nChInFile = props.numChan; %total number of recorded channels
dataType = 'int16';

%%%%% create memory map %%%%%
mmf = memmapfile([kilosort_path '\' dir([kilosort_path '\' sprintf('Neuropixels_%s*',subj)]).name], 'Format', {dataType, [nChInFile nSamp], 'x'});

%%%%% loop through clusters to create metrics struct %%%%%
metrics = [];
for clu = 1:length(rawclusters{1}) %each session should have same number of clusters, and already only good clusters
    %get waveforms from raw data (from makeWFstructure_fromRaw)
    maxChan = rawclusters{1}(clu).maxChan + 1; %max amplitude channel (0-based)
    tmp_nWFsToLoad = min(nWFsToLoad, length(rawclusters_allrec(clu).spikeInds));
    % tmp_nWFsToLoad = length(theseST); %load all spikes detected from the unit
    extractST = round(rawclusters_allrec(clu).spikeInds(randperm(length(rawclusters_allrec(clu).spikeInds), tmp_nWFsToLoad)));
    nSamps2Chk = round(0.2 / (1/props.sampRate *1000));

    theseWF = nan(tmp_nWFsToLoad, nWFsamps);
    for ii=1:tmp_nWFsToLoad
        %eliminate spikes less than 1s into rec file && index  + 2ms after
        %spike does not go over the data boundary
        if (extractST(ii)/props.sampRate > 1.0) && ((extractST(ii) + tAroundSpike(2)) < nSamp)

            %look for minimum amplitude +/- 0.2ms around the K2 spikeIdx
            a = extractST(ii) - nSamps2Chk : extractST(ii) + nSamps2Chk;
            tempWF = mmf.Data.x(maxChan, a);

            b = find(tempWF == min(tempWF));
            minIdx = a(b); %index of minimum waveform amplitude

            %only load from the max amplitude channel
            theseWF(ii,:) = mmf.Data.x(maxChan, minIdx-tAroundSpike(1):minIdx+tAroundSpike(2));
        end
    end
    %subtract from baseline to go from zero to peak spike ampltiude as trough
    theseWF = theseWF - theseWF(:,1);

    %get isi info
    isi = []; spikeCount = [];
    for r = 1:size(props.fileNames,1)
        isi = [isi; diff(rawclusters{r}(clu).spikeInds)];
        spikeCount = [spikeCount; length(rawclusters{r}(clu).spikeInds)];
    end

    %create metrics struct
    spikeCount = sum(spikeCount);
    recLengthAll = nSamp./props.sampRate; %in s
    metrics(clu).ID = rawclusters{1}(clu).ID;
    metrics(clu).WF.mn = mean(theseWF, 1, 'omitnan');
    metrics(clu).WF.std = std(theseWF, 0, 1, 'omitnan');
    metrics(clu).snr = (max(metrics(clu).WF.mn) - min(metrics(clu).WF.mn))/mean(metrics(clu).WF.std);
    metrics(clu).WF.info = '1ms before and 2ms after minimum amplitude';
    metrics(clu).firingrate = spikeCount/recLengthAll; %in Hz
    metrics(clu).isi.all_ms = isi./(props.sampRate/1000); %in ms
    metrics(clu).isi.edges_ms = 0:0.1:10;
    metrics(clu).isi_h = nan(1, length(metrics(clu).isi.edges_ms));
    metrics(clu).isi_h(1:length(metrics(clu).isi.edges_ms)) = histc(metrics(clu).isi.all_ms, 0:0.1:10);
    metrics(clu).numspikes = spikeCount;
    metrics(clu).files = props.fileNums;
    metrics(clu).index = [subj '_' sessDate];
    metrics(clu).maxChan = rawclusters{1}(clu).maxChan;
    metrics(clu).samprate = props.sampRate;
    %(end makeWFstructure_fromRaw)

    %find stable clusters and add info to metrics struct (from getstableclustertimes_gauss_K2)
    plotexamples = 0;
    totaltime = 0;
    mintimestable = 5;

    %get spikes and recording file durations
    for r = 1:size(props.fileNames,1)
        if r == 1
            fr = [];
            fr.totalspiketimes = [];
            fr.totalspiketimes = rawclusters{r}(clu).spikeInds'./(props.sampRate/1000); %put into ms
        else
            fr.totalspiketimes = [fr.totalspiketimes (rawclusters{r}(clu).spikeInds'./(props.sampRate/1000)+totaltime)]; %all spike times
        end

        totaltime = totaltime + props.recLength(r)/props.sampRate*1000; %total time of the first recording

        if r == 1
            switchdur = [0 props.recLength(r)/props.sampRate*1000];
        else
            switchdur = [switchdur; (switchdur(end,2)+1) (switchdur(end,2)+1+props.recLength(r)/props.sampRate*1000)]; %total time of recordings after the first
        end
    end

    %gauss FR
    spiketrainedges = 0:10:totaltime; %10 ms bins
    fr.spiketrain = histc(fr.totalspiketimes, spiketrainedges);
    fr.gaussw = gausswin(10*100*60)*100./(sum(gausswin(10*100*60))); %area is 100 to get to Hz
    fr.gaussfr = conv(fr.spiketrain, fr.gaussw, 'same');

    groups = kmeans(fr.gaussfr', 2, 'MaxIter', 1000); %assume FR bimodal, get kmeans %changed MaxIter to 1000 from default ALP 11/28

    [amean, id] = min([mean(fr.gaussfr(groups==1)), mean(fr.gaussfr(groups==2))]); %low mean
    bmean = max([mean(fr.gaussfr(groups==1)), mean(fr.gaussfr(groups==2))]); %high mean
    astd = std(fr.gaussfr(groups==id));

    %apply kmeans thresholding if the FR drops to 10% of high FR mean
    if min(fr.gaussfr(100*60:end-100*60)) < 0.1*bmean %to try and help with edge effects
        findtimes = find(fr.gaussfr >= 2*astd+amean); %threshold is 2std above low mean
        if length(findtimes) > 1 && sum(diff(findtimes) == 1) > 0
            temptimes = contiguous(diff(findtimes),1);
            temptimes = temptimes{2};
            tempdiff = temptimes(:,2) - temptimes(:,1);
            binid = find(tempdiff == max(tempdiff)); %max duration above threshold

            %GaussFR
            fr.incltimes = [spiketrainedges(findtimes(temptimes(binid,1))) spiketrainedges(findtimes(temptimes(binid,2)+1))]; %in ms
            findfilesmat = [isExcluded(switchdur(:,1),fr.incltimes) isExcluded(switchdur(:,2),fr.incltimes)];

            if (fr.incltimes(2) - fr.incltimes(1))/(1000*60) >= mintimestable
                for file = 1:size(findfilesmat,1)
                    if sum(findfilesmat(file,:)) > 0
                        fr.stabletimes(file,:) = [max(fr.incltimes(1), switchdur(file,1)) min(fr.incltimes(2), switchdur(file,2))];
                        fr.stabletimes(file,:) = fr.stabletimes(file,:) - switchdur(file,1)*ones(1,2);
                        %added ALP 4/10/18
                    elseif sum(findfilesmat(:,1)) == 0 && sum(findfilesmat(:,2)) == 0
                        temponerec_1 = isExcluded(fr.incltimes(1), switchdur(file,:));
                        temponerec_2 = isExcluded(fr.incltimes(2), switchdur(file,:));

                        %if the stable time only falls into one recording
                        % and if that recording is this file, then make
                        % stable times, otherwise set to 0 for this file
                        if isequal(temponerec_1, temponerec_2) && sum(temponerec_1+temponerec_2) > 0
                            fr.stabletimes(file,:) = [max(fr.incltimes(1), switchdur(file,1)) min(fr.incltimes(2), switchdur(file,2))];
                            fr.stabletimes(file,:) = fr.stabletimes(file,:) - switchdur(file,1)*ones(1,2);
                        else
                            fr.stabletimes(file,:) = [0 0];
                        end
                        %ALP 4/10 end
                    else
                        fr.stabletimes(file,:) = [0 0];
                    end
                end
            else
                fr.stabletimes = zeros(size(findfilesmat,1),2);
                fr.incltimes = [0 0];
            end
        else
            fr.stabletimes = zeros(size(findfilesmat,1),2);
            fr.incltimes = [0 0];
        end

        if plotexamples
            plottime = spiketrainedges./1000/60;
            %patch for plotting
            patchindsx = [fr.incltimes(1)./(1000*60) fr.incltimes(1)/(1000*60) fr.incltimes( 2)/(1000*60) fr.incltimes(2)/(1000*60)];
            patchindsy = [0 (max(fr.gaussfr)+0.1*max(fr.gaussfr)) (max(fr.gaussfr)+0.1*max(fr.gaussfr)) 0];

            figure
            patch(patchindsx, patchindsy, [0.6 0.8 1.0])
            hold on
            alpha(0.1)
            ylim([0 (max(fr.gaussfr)+0.1*max(fr.gaussfr))])
            plot(plottime, fr.gaussfr, 'k', 'LineWidth', 1.5)
            hold on
            plot(plottime, (amean+2*astd)*ones(1,length(fr.gaussfr)), 'r--', 'LineWidth', 1.5)
            plot(plottime, 0.1*bmean*ones(1,length(fr.gaussfr)), 'b--', 'LineWidth', 1.5)
            legend('Included', 'FR', '2std', '10% High')
            ylabel('Firing Rate (Hz)')
            xlabel('Time (min)')
            xlim([0 totaltime/(1000*60)])
        end

    else %else whole time stable
        for file = 1:size(props.fileNames,1)
            fr.incltimes = [0 totaltime];
            fr.stabletimes(file,:) = [0 (switchdur(file,2)-switchdur(file,1))];
        end

        if plotexamples
            plottime = spiketrainedges./1000/60;

            figure
            hold on
            plot(plottime, fr.gaussfr, 'k')
            hold on
            plot(plottime, 0.1*bmean*ones(1,length(fr.gaussfr)), 'b--')
            legend('Firing Rate', '10% High FR')
            ylabel('FR (Hz)')
            xlabel('Time (min)')
            xlim([0 totaltime/(1000*60)])
        end
    end

    fr.stabletimes = fr.stabletimes./1000; %ms to s
    stabletimes = fr.stabletimes;

    %get mean and peak FR
    stableidx = fr.incltimes/10; %fr.incltimes in ms, fr.gausswin in ms/10 (10ms bins)

    %round to integers to get indices
    stableidx = round(stableidx);
    if stableidx(1) == 0
        stableidx(1) = 1;
    end

    %calc mean/peak FR
    if stableidx(2) ~= 0 %accounts for cluster with no stable times
        stableFR = fr.gaussfr(stableidx(1):stableidx(2));
        peakFR = max(stableFR);
        meanFR= mean(stableFR);
    else
        peakFR = nan;
        meanFR = nan;
    end

    %save to metrics struct
    metrics(clu).stable.times = stabletimes;
    metrics(clu).stable.meanFR = meanFR;
    metrics(clu).stable.peakFR = peakFR;
    metrics(clu).stable.info = 'times in [s], meanFR and peakFR from stable period over all recordings, in [Hz], idx in samprate';
    %(end getstableclustertimes_gauss_K2)

end%clu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% quality control to find good clusters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Quality control thresholds (from Kilosort2_Pipeline_Abby.m) %%%%%
th.SNR =  1;                    % >= 1 SNR
th.ISI = 0.008;                 % <= 0.8% refractory period violations
th.refractoryPeriod = 0.001;    % 1ms refractory period duration
th.info = '>= th.SNR, <= th.ISI (frac violations/allISI), th.refractoryPeriod in s';

%%%%% get good cells %%%%%
%adapted from applyQualityMetrics.m

%SNR
good_snr = [metrics.snr] >= th.SNR;
snrexcl = sum(double(good_snr) == 0);
disp([num2str(snrexcl), ' of ', num2str(length(metrics)), ' excluded for SNR < ', num2str(th.SNR)])

%refractory period violations
isibins = metrics(1).isi.edges_ms < th.refractoryPeriod*1000; %can use the first ind bc all edges the same
temp = reshape([metrics.isi_h], [length(isibins), length(metrics)])'; %get matrix of isi hist of all clusters
temp = sum(temp(:,isibins),2); %get the number of spikes in the refractory period for each cluster
temp = temp'./[metrics.numspikes]; %normalize by total spike count for each cluster
good_isi = temp <= th.ISI;
isiexcl = sum(double(good_isi) == 0);
disp([num2str(isiexcl), ' of ', num2str(length(metrics)), ' excluded for > ' num2str(th.ISI*100),'% refractory violations'])

%combine all metrics
temp = zeros(1, length(metrics));
temp((good_snr & good_isi)) = 1; %fixed bug when two logical values are equally zero, it gets included as a good unit when it shouldn't: NJ 19.09.06
good_final = temp; %vector of 0 and 1s, 1x(#cells), 0 = excluded cell, 1 = included cell

%display information about how many cells survived the quality metrics into
%the command line
totalexcl = sum(double(good_final) == 0);
totalincl = sum(double(good_final) == 1);
disp([num2str(totalexcl), ' of ', num2str(length(metrics)), ' clusters excluded.'])
disp([num2str(totalincl), ' clusters survived.'])

%%%%% make good clusters structure %%%%%
clusters = [];
for r = 1:size(props.fileNames,1)
    %populate structure
    clusters(r).samprate = props.sampRate;
    clusters(r).th = th;
    clusters(r).info = 'single units in .data struct. spikeInds is spike index with samprate denoted by field .samprate. maxChan is 0 based.';
    fileNum = props.fileNums(r);
    clusters(r).index = [subj '_' sessDate '_' num2str(fileNum)];
    numgood = 1;

    for clu = 1:length(good_final)
        if good_final(clu)
            clusters(r).data(numgood).ID = rawclusters{r}(clu).ID;
            clusters(r).data(numgood).maxChan = rawclusters{r}(clu).maxChan;
            clusters(r).data(numgood).spikeInds = rawclusters{r}(clu).spikeInds;
            clusters(r).data(numgood).spikeAmps = rawclusters{r}(clu).spikeAmps;
            numgood = numgood+1;
        end
    end
end

%%%%% make good cluster metrics structure %%%%%
clustermetrics = metrics(logical(good_final));

%%%%% make good clusters structure for all recordings %%%%%
clusters_allrec = rawclusters_allrec(logical(good_final));
[clusters_allrec.info] = deal({'all files. post quality control metrics'});
%(end applyQualityMetrics.m)

%%%%% ensure 100% every cluster is "good" %%%%%
%Note: JLK 8/5/25 having issues with some clusters included in clusters.mat
%not actually good clusters, leading to errors later when compairing clusters and
%cell_metrics. These lines 100% ensures only good clusters are saved.
for r = 1:size(props.fileNames,1)
    clusters(r).data = clusters(r).data(ismember([clusters(r).data.ID], goodUnits));
end
clustermetrics = clustermetrics(ismember([clusters(1).data.ID], goodUnits));
clusters_allrec = clusters_allrec(ismember([clusters(1).data.ID], goodUnits));

%save the good clusters + metrics to kilosort folder
save([neuralRawDataPath '\kilosort4\clusters.mat'], 'clusters','-v7.3')
save([neuralRawDataPath '\kilosort4\clustermetrics.mat'], 'clustermetrics','-v7.3')
save([neuralRawDataPath '\kilosort4\clusters_allrec.mat'], 'clusters_allrec','-v7.3')

end%function