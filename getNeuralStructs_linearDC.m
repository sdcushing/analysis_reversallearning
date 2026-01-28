function [rawDataBySessionNeural, rawDataByLapNeural, rawDataByTrialNeural] = getNeuralStructs_linearDC(subj,sessDate,sessNum,params,virmenSessDataPath, processedDataPath, saveNeuralPath)
%adapted from getRestSessionStats_linearJLK
%
%% Create Raw Data Struct for Whole Session with LFP and AP Times and Data %%
%need ttls , virmen , lfp, spikes (already clustered)
%lfp is raweeg in processed data (ind. folders for regions + channels)
%[processed data path + ]
%spikes from kilosort/cellexplorer [processeddatapath]
%virmen from virmen (databysession in output structs) [virmensessdatapath]
%ttl from digital input (nidaq for np, reward file in processed data for
%intan) [processeddatapath]

%%i think best plan is to use rawpos, which already has ephys ind all
%%aligned, instead of JK's behavior structure. since it's ephysind, not ts,
%%will need to make sure I have appropriate indexes. need to look through.
%%raweeg + velocity + has same samples as sync, and clusters has same
%%samprate, so actually should be all good.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load raw virmen data %%%%%switch to rawpos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawDataBySession = load([virmenSessDataPath '\rawDataBySession.mat']);%whole session
rawDataBySession = rawDataBySession.rawDataBySession;
%%%%% deal with special cases %%%%%leaving these in so I have ideas how to
%%%%% deal with them
%no start TTL
% Assumption: if first TTL is past halfway through session, assume missed start TTL
% Solution: subtract 20 minutes from last TTL
% if virmen_ttl_start_ind > length(rawDataBySession.vrTime)/2 
%     %subtract 20 minutes from last TTL and new start index and time
%     %NIDAQ
%     nidaq_ttl_start_ind = nidaq_ttl_end_ind - (20 * 60 * nidaq_samprate);
%     nidaq_ttl_start_time = nidaq_ttl_start_ind/nidaq_samprate;%seconds
%     nidaq_ttl_duration = nidaq_ttl_end_time - nidaq_ttl_start_time;%seconds
% 
%     %VIRMEN
%     tmp_virmen_start_time = virmen_ttl_end_time - (20 * 60);%seconds
%     tmp_virmen_rec_time_diffs = [0; diff(rawDataBySession.vrTime)];
%     tmp_virmen_rec_time_sums = cumsum(tmp_virmen_rec_time_diffs);
%     virmen_ttl_start_ind = find(tmp_virmen_rec_time_sums>=tmp_virmen_start_time,1,'first');
%     virmen_ttl_start_time = sum([0; diff(rawDataBySession.vrTime(1:virmen_ttl_start_ind))]);%seconds
%     virmen_ttl_duration = virmen_ttl_end_time - virmen_ttl_start_time;%seconds
% 
%     %Duration check
%     if abs(virmen_ttl_duration - nidaq_ttl_duration) > 1%seconds
%         sprintf('Check duration of TTL pulses for Virmen vs. NIDAQ.')
%     else
%         sprintf('Resolved no start TTL pulse.')
%     end%duration check
% end%no start TTL
% 
% %dropped data
% % Solution (non-optimal): cut data before error message time
% if length(find(diff(nidaq_data(1,:))>4)) ~= length(find(rawDataBySession.ttl==1))
%     if strcmp(subj,'JK15') && strcmp(sessDate,'250530') && strcmp(sessNum,'3')
%         %add 15 minutes from start TTL and end start index and time
%         %NIDAQ
%         nidaq_ttl_end_ind = nidaq_ttl_start_ind + (15 * 60 * nidaq_samprate);
%         nidaq_ttl_end_time = nidaq_ttl_end_ind/nidaq_samprate;%seconds
%         nidaq_ttl_duration = nidaq_ttl_end_time - nidaq_ttl_start_time;%seconds
% 
%         %VIRMEN
%         tmp_virmen_end_time = virmen_ttl_start_time + (15 * 60);%seconds
%         tmp_virmen_rec_time_diffs = [0; diff(rawDataBySession.vrTime)];
%         tmp_virmen_rec_time_sums = cumsum(tmp_virmen_rec_time_diffs);
%         virmen_ttl_end_ind = find(tmp_virmen_rec_time_sums<=tmp_virmen_end_time,1,'last');
%         virmen_ttl_end_time = sum([0; diff(rawDataBySession.vrTime(1:virmen_ttl_end_ind))]);%seconds
%         virmen_ttl_duration = virmen_ttl_end_time - virmen_ttl_start_time;%seconds
% 
%         %Duration check
%         if abs(virmen_ttl_duration - nidaq_ttl_duration) > 1%seconds
%             sprintf('Check duration of TTL pulses for Virmen vs. NIDAQ.')
%         else
%             sprintf('Resolved discrepency in Virmen vs. NIDAQ TTL pulses.')
%         end%duration check
%     end%JK15_250530_3
% end%different number of ttl pulses

%%%%% LFPs %%%%%

lfp_filepath = fullfile(string(processedDataPath), string(params.brainReg));
% Get list of subfolders
d = dir(lfp_filepath);
isSubfolder = [d.isdir];
names = {d(isSubfolder).name};
names = names(~ismember(names, {'.','..'}));   % remove . and ..
numFolders = numel(names);

% Preallocate after loading first file to know sample count
firstFile = fullfile(lfp_filepath, names{1}, sprintf('raweeg%s.mat', sessNum));
tmp = load(firstFile);
samples = numel(tmp.raweeg.data);
lfp_data = zeros(numFolders, samples/10);%want size to fit downsampled data (20k/2k)
% Initialize metadata structure
lfp_meta = struct();
lfp_meta.samprate = tmp.raweeg.samprate;   % same for all files
% Preallocate per‑channel metadata fields
lfp_meta.port          = cell(numFolders,1);
lfp_meta.nativechannel = cell(numFolders,1);
lfp_meta.channelInd    = zeros(numFolders,1);
lfp_meta.gapInfo       = cell(numFolders,1);
lfp_samprate_down = 2000;

% Loop through folders
for i = 1:numFolders
    folderPath = fullfile(lfp_filepath, names{i});
    matFile = fullfile(folderPath, sprintf('raweeg%s.mat', sessNum));

    S = load(matFile);
    R = S.raweeg;

    %filter EEG
    [eeg, outlierthresh, numoutlierperiods, outlierindices, outlierperiods] = makefiltereeg_Intan_DC(R, params.filteegfreq, params.eegsamprate, params.outliernstd);

    % Store LFP data
    lfp_data(i, :) = double(eeg.data);

    % Store metadata
    lfp_meta.port{i}          = R.port;
    lfp_meta.nativechannel{i} = R.nativechannel;
    lfp_meta.channelInd(i)    = R.channelInd;
    lfp_meta.gapInfo{i}       = R.gapInfo;
    lfp_meta.outlierthresh{i}    = outlierthresh;
    lfp_meta.numoutlierperiods{i} = numoutlierperiods;
    lfp_meta.outlierindices{i} = outlierindices;
    lfp_meta.outlierperiods{i} = outlierperiods;
    
end
%will want to save lfp.meta somewhere as well
lfp_meta.downsamplerate = lfp_samprate_down;
%downsample from 20000 to 2000 Hz for ripple filter later

% lfp_data = lfp_data(:,1:downSamp:end);

%need to filter outliers + downsample instead. filter works on individual
%channels, need to rework filter. makefiltereeg_Intan,
%interpoveroutliers_041819, filtereeg


%remove 60 Hz noise; optional because takes long time and lfp trace from example session
% looked similar before and after
if params.rm60HzNoise
    lfp_data_no60 = [];
    chparams.fpass=[1 400]; % band of frequencies to be kept
    chparams.tapers=[6 11]; % taper parameters
    chparams.pad=1; % pad factor for fft
    chparams.trialave=1;
    chparams.Fs = lfp_samprate_down;
    for ch = 1:size(lfp_data,1)
        lfp_data_no60(ch,:) = rmlinesmovingwinc(lfp_data(ch,:),[4 2],10,chparams,.00000001,'n', 60);
    end%ch
    lfp_data = lfp_data_no60;
end%params.rm60HzNoise

%determine lfp start/end times
virmen_rec_time_lfp = round(rawDataBySession.ephysInd/lfp_samprate_down);%20k/2k = 10 (og samprate/downsamp)

% %ensure last virmen data index is around last lfp ttl index
% if abs( (lfp_ttl_end_ind/lfp_samprate_down) - (virmen_rec_time_lfp(end)/lfp_samprate_down) ) > 1%seconds
%     sprintf('Check conversion of virmen times to lfp times.')
% end%time check


%%%%% APs %%%%%
%load session clusters struct
load([processedDataPath '\kilosort4\clusters.mat'])
load([processedDataPath '\kilosort4\' sprintf('%s_%s.cell_metrics.cellinfo.mat', subj, sessDate)])

%extract ap data
%clustersRow = []; for cl = 1:length(clusters); clustersRow(cl) = strcmp(clusters(cl).index,sprintf('%s_%s_%s',subj,sessDate,sessNum));end%cl
ap_samprate = clusters.samprate;
ap_data = clusters.data;

virmen_rec_time_ap = rawDataBySession.ephysInd;

% %ensure last virmen data index is around last ap ttl index
% if abs( (ap_ttl_end_ind/ap_samprate) - (virmen_rec_time_ap(end)/ap_samprate) ) > 1%seconds
%     sprintf('Check conversion of virmen times to lfp times.')
% end%time check

%remove spikes outside of virmen session
ap_ttl_start_ind = virmen_rec_time_ap(1); ap_ttl_end_ind = virmen_rec_time_ap(end);
for clu = 1:length(ap_data)
    ap_data(clu).spikeInds = ap_data(clu).spikeInds(ap_data(clu).spikeInds >= ap_ttl_start_ind & ap_data(clu).spikeInds <= ap_ttl_end_ind);
    ap_data(clu).spikeAmps = ap_data(clu).spikeAmps(ap_data(clu).spikeInds >= ap_ttl_start_ind & ap_data(clu).spikeInds <= ap_ttl_end_ind);
end%clu

%add putative cell type
for clu = 1:length(ap_data)
    ap_data(clu).putativeCellType = cell_metrics.putativeCellType{cell_metrics.cluID == ap_data(clu).ID};    
end%clu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% create sturct and save data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawDataBySessionNeural = [];
rawDataBySessionNeural.vrTime = rawDataBySession.vrTime;
rawDataBySessionNeural.lfpTime = virmen_rec_time_lfp;
rawDataBySessionNeural.apTime = virmen_rec_time_ap;
rawDataBySessionNeural.lfpData = lfp_data;%already reduced to start to end TTL for finding outliers above
rawDataBySessionNeural.lfpmeta = lfp_meta;
% rawDataBySessionNeural.lfpOutlierInd = [tmpHistStInds'-2 tmpHistStInds'+7];
rawDataBySessionNeural.apData = ap_data;
rawDataBySessionNeural.currentDeg = rawDataBySession.currentDeg;
rawDataBySessionNeural.speed = [0; diff(rawDataBySessionNeural.currentDeg) ./ diff(rawDataBySessionNeural.vrTime)];
rawDataBySessionNeural.speedSmooth = gaussSmooth(rawDataBySessionNeural.speed', 5)';
rawDataBySessionNeural.isMoving = rawDataBySessionNeural.speed > params.speedTh;
rawDataBySessionNeural.rewarded = [0; diff(rawDataBySession.rewards)];
rawDataBySessionNeural.licked = [0; diff(rawDataBySession.licks)];
rawDataBySessionNeural.currentZone = rawDataBySession.currentZone;

filename = [saveNeuralPath '\' 'rawDataBySessionNeural.mat'];
save(filename, 'rawDataBySessionNeural', '-v7.3');

%% Create Raw Data Struct for Laps with LFP and AP Times and Data %%
%Note: Rest sessions do not have lap structs

if isfile([virmenSessDataPath '\rawDataByLap.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% load raw virmen data %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rawDataByLap = load([virmenSessDataPath '\rawDataByLap.mat']);%laps
    rawDataByLap = rawDataByLap.rawDataByLap;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% create sturct %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    rawDataByLapNeural = [];
    for ii = 1:length(rawDataByLap)%loop through laps

        %%%%% index whole session neural struct using virmen lap times %%%%%
        virmen_lap_start_ind = find(rawDataBySessionNeural.vrTime >= rawDataByLap(ii).vrTime(1),1,'first');
        virmen_lap_end_ind = find(rawDataBySessionNeural.vrTime <= rawDataByLap(ii).vrTime(end),1,'last');

        %%%%% add data to struct %%%%%
        rawDataByLapNeural(ii).vrTime = rawDataBySessionNeural.vrTime(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).lfpTime = rawDataBySessionNeural.lfpTime(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).apTime = rawDataBySessionNeural.apTime(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).lfpData = rawDataBySessionNeural.lfpData(:,(rawDataByLapNeural(ii).lfpTime(1)):(rawDataByLapNeural(ii).lfpTime(end)));
        for clu = 1:length(rawDataBySessionNeural.apData)
            rawDataByLapNeural(ii).apData(clu).ID = rawDataBySessionNeural.apData(clu).ID;
            rawDataByLapNeural(ii).apData(clu).maxChan = rawDataBySessionNeural.apData(clu).maxChan;
            rawDataByLapNeural(ii).apData(clu).spikeInds = rawDataBySessionNeural.apData(clu).spikeInds(rawDataBySessionNeural.apData(clu).spikeInds >= rawDataByLapNeural(ii).apTime(1) & rawDataBySessionNeural.apData(clu).spikeInds <= rawDataByLapNeural(ii).apTime(end));
            rawDataByLapNeural(ii).apData(clu).spikeAmps = rawDataBySessionNeural.apData(clu).spikeAmps(rawDataBySessionNeural.apData(clu).spikeInds >= rawDataByLapNeural(ii).apTime(1) & rawDataBySessionNeural.apData(clu).spikeInds <= rawDataByLapNeural(ii).apTime(end));
            rawDataByLapNeural(ii).apData(clu).putativeCellType = rawDataBySessionNeural.apData(clu).putativeCellType;
        end%clu
        rawDataByLapNeural(ii).currentDeg = rawDataBySessionNeural.currentDeg(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).speed = rawDataBySessionNeural.speed(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).speedSmooth = rawDataBySessionNeural.speedSmooth(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).isMoving = rawDataBySessionNeural.isMoving(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).rewarded = rawDataBySessionNeural.rewarded(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).licked = rawDataBySessionNeural.licked(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).currentZone = rawDataBySessionNeural.currentZone(virmen_lap_start_ind:virmen_lap_end_ind);

        %%%%% binned data for decoding %%%%%
        %DISTANCE%
        %create bins and identify which bins each position data point belongs to
        degBinEdges = 0:params.binsize_deg:360;
        rawDataByLapNeural(ii).degBinEdges = degBinEdges;
        degBinIden = discretize(rawDataByLapNeural(ii).currentDeg,degBinEdges);
        rawDataByLapNeural(ii).degBinSize = params.binsize_deg;

        %compute occupancy, speed, and licks per spatial bin
        for b = 1:length(degBinEdges)-1
            rawDataByLapNeural(ii).degBinOccup(b) = sum(diff(rawDataByLapNeural(ii).vrTime(degBinIden == b)));
            rawDataByLapNeural(ii).degBinPos(b) = sum(diff(rawDataByLapNeural(ii).currentDeg(degBinIden == b)));
            rawDataByLapNeural(ii).degBinLicked(b) = sum(rawDataByLapNeural(ii).licked(degBinIden == b));
        end%bin
        rawDataByLapNeural(ii).degBinSpeed = rawDataByLapNeural(ii).degBinPos ./ rawDataByLapNeural(ii).degBinOccup;
        rawDataByLapNeural(ii).degBinLickRate = rawDataByLapNeural(ii).degBinLicked ./ rawDataByLapNeural(ii).degBinOccup;
        %smooth data
        rawDataByLapNeural(ii).degBinOccupSmooth = gaussSmooth(rawDataByLapNeural(ii).degBinOccup, 2)';
        rawDataByLapNeural(ii).degBinPosSmooth = gaussSmooth(rawDataByLapNeural(ii).degBinPos, 2)';
        rawDataByLapNeural(ii).degBinSpeedSmooth = gaussSmooth(rawDataByLapNeural(ii).degBinSpeed, 2)';
        rawDataByLapNeural(ii).degBinLickRateSmooth = gaussSmooth(rawDataByLapNeural(ii).degBinLickRate, 2)';

        %add spiking data and rate map
        rawDataByLapNeural(ii).degBinSpikeCount = zeros(length(degBinEdges)-1,length(rawDataBySessionNeural.apData));
        rawDataByLapNeural(ii).degBinSpikeCountSmooth = zeros(length(degBinEdges)-1,length(rawDataBySessionNeural.apData));
        rawDataByLapNeural(ii).degBinRateMap = zeros(length(degBinEdges)-1,length(rawDataBySessionNeural.apData));
        for clu = 1:length(rawDataBySessionNeural.apData)
            spikePosInd = lookup2(rawDataByLapNeural(ii).apData(clu).spikeInds, rawDataByLapNeural(ii).apTime);
            spikePos = rawDataByLapNeural(ii).currentDeg(spikePosInd);
            rawDataByLapNeural(ii).degBinSpikeCount(:,clu) = histcounts(spikePos,degBinEdges);
            rawDataByLapNeural(ii).degBinSpikeCountSmooth(:,clu) = gaussSmooth(rawDataByLapNeural(ii).degBinSpikeCount(:,clu)', 2)';
            rawDataByLapNeural(ii).degBinRateMap(:,clu) = rawDataByLapNeural(ii).degBinSpikeCountSmooth(:,clu) ./ rawDataByLapNeural(ii).degBinOccupSmooth;
        end%clu

    end%ii

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% save data %%%%%
    %%%%%%%%%%%%%%%%%%%%%
    filename = [saveNeuralPath '\' 'rawDataByLapNeural.mat'];
    save(filename, 'rawDataByLapNeural', '-v7.3');

end%if isfile([virmenSessDataPath '\rawDataByLap.mat'])


%% Create Raw Data Struct for Trials with LFP and AP Times and Data %%
%Note: Rest sessions and active sessions with <= 1 trial do not have trial structs

if isfile([virmenSessDataPath '\rawDataByTrial.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% load raw virmen data %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rawDataByTrial = load([virmenSessDataPath '\rawDataByTrial.mat']);%trials
    rawDataByTrial = rawDataByTrial.rawDataByTrial;

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% create sturct %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    rawDataByTrialNeural = cell(size(rawDataByTrial));
    %loop through trials for each zone
    for znType = 1:size(rawDataByTrial,1) %1=reward, 2=nonreward, 3 = alt nonreward
        for znNum = 1:size(rawDataByTrial,2)-1%always 3 types of zones
            for tr = 1:size(rawDataByTrial{znType,znNum},2)

                %%%%% index whole session neural struct using virmen trial times %%%%%
                if isempty(rawDataByTrial{znType,znNum})
                    continue
                end                
                virmen_trial_start_ind = find(rawDataBySessionNeural.vrTime >= rawDataByTrial{znType,znNum}(tr).vrTime(1),1,'first');
                virmen_trial_end_ind = find(rawDataBySessionNeural.vrTime <= rawDataByTrial{znType,znNum}(tr).vrTime(end),1,'last');

                %%%%% add data to struct %%%%%
                rawDataByTrialNeural{znType,znNum}(tr).vrTime = rawDataBySessionNeural.vrTime(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).lfpTime = rawDataBySessionNeural.lfpTime(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).apTime = rawDataBySessionNeural.apTime(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).lfpData = rawDataBySessionNeural.lfpData(:,(rawDataByTrialNeural{znType,znNum}(tr).lfpTime(1)):(rawDataByTrialNeural{znType,znNum}(tr).lfpTime(end)));
                for clu = 1:length(rawDataBySessionNeural.apData)
                    rawDataByTrialNeural{znType,znNum}(tr).apData(clu).ID = rawDataBySessionNeural.apData(clu).ID;
                    rawDataByTrialNeural{znType,znNum}(tr).apData(clu).maxChan = rawDataBySessionNeural.apData(clu).maxChan;
                    rawDataByTrialNeural{znType,znNum}(tr).apData(clu).spikeInds = rawDataBySessionNeural.apData(clu).spikeInds(rawDataBySessionNeural.apData(clu).spikeInds >= rawDataByTrialNeural{znType,znNum}(tr).apTime(1) & rawDataBySessionNeural.apData(clu).spikeInds <= rawDataByTrialNeural{znType,znNum}(tr).apTime(end));
                    rawDataByTrialNeural{znType,znNum}(tr).apData(clu).putativeCellType = rawDataBySessionNeural.apData(clu).putativeCellType;
                end%clu
                rawDataByTrialNeural{znType,znNum}(tr).currentDeg = rawDataBySessionNeural.currentDeg(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).speed = rawDataBySessionNeural.speed(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).speedSmooth = rawDataBySessionNeural.speedSmooth(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).rewarded = rawDataBySessionNeural.rewarded(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).isMoving = rawDataBySessionNeural.isMoving(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).licked = rawDataBySessionNeural.licked(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).currentZone = rawDataBySessionNeural.currentZone(virmen_trial_start_ind:virmen_trial_end_ind);

                %%%%% binned data for decoding %%%%%
                %DISTANCE%
                %create bins and identify which bins each position data point belongs to
                if floor(rawDataByTrialNeural{znType,znNum}(tr).currentDeg(1)) + (params.gapAfter*params.binsize_deg + params.cueSize) > 360 %last bin > 360 deg, use wrapTo360
                    degBinEdges = wrapTo360(floor(rawDataByTrialNeural{znType,znNum}(tr).currentDeg(1)):params.binsize_deg:ceil(rawDataByTrialNeural{znType,znNum}(tr).currentDeg(end)+360));
                    if isempty(degBinEdges)
                        error('degBinEdges is empty, stopping execution.');
                    end
                    rawDataByTrialNeural{znType,znNum}(tr).degBinEdges = degBinEdges; %save to struct before adding 360 to make edges monotonically increasing for discretize function
                    transition = find(diff(degBinEdges)<-300) + 1;%find index for the start of the next lap
                    degBinEdges(transition:end) = degBinEdges(transition:end)+360;%add 360 to make edges monotomically increase for discretize function
                    tmpCurrentDeg = rawDataByTrialNeural{znType,znNum}(tr).currentDeg;
                    transition = find(diff(tmpCurrentDeg)<-300) + 1;
                    tmpCurrentDeg(transition:end) = rawDataByTrialNeural{znType,znNum}(tr).currentDeg(transition:end)+360;
                    try
                        degBinIden = discretize(tmpCurrentDeg,degBinEdges);
                    catch
                        error('Issue discretizing current degrees into degree bins')
                    end
                else%last bin <= 360 deg
                    degBinEdges = floor(rawDataByTrialNeural{znType,znNum}(tr).currentDeg(1)):params.binsize_deg:ceil(rawDataByTrialNeural{znType,znNum}(tr).currentDeg(end));
                    if isempty(degBinEdges)
                        error('degBinEdges is empty, stopping execution.');
                    end
                    rawDataByTrialNeural{znType,znNum}(tr).degBinEdges = degBinEdges;
                    degBinIden = discretize(rawDataByTrialNeural{znType,znNum}(tr).currentDeg,degBinEdges);
                end
                rawDataByTrialNeural{znType,znNum}(tr).degBinSize = params.binsize_deg;

                %compute occupancy, speed, and licks per spatial bin 
                for b = 1:length(degBinEdges)-1
                    rawDataByTrialNeural{znType,znNum}(tr).degBinOccup(b) = sum(diff(rawDataByTrialNeural{znType,znNum}(tr).vrTime(degBinIden == b)));
                    rawDataByTrialNeural{znType,znNum}(tr).degBinPos(b) = sum(diff(rawDataByTrialNeural{znType,znNum}(tr).currentDeg(degBinIden == b)));
                    rawDataByTrialNeural{znType,znNum}(tr).degBinLicked(b) = sum(rawDataByTrialNeural{znType,znNum}(tr).licked(degBinIden == b));
                end%bin
                rawDataByTrialNeural{znType,znNum}(tr).degBinSpeed = rawDataByTrialNeural{znType,znNum}(tr).degBinPos ./ rawDataByTrialNeural{znType,znNum}(tr).degBinOccup;
                rawDataByTrialNeural{znType,znNum}(tr).degBinLickRate = rawDataByTrialNeural{znType,znNum}(tr).degBinLicked ./ rawDataByTrialNeural{znType,znNum}(tr).degBinOccup;
                %smooth data
                rawDataByTrialNeural{znType,znNum}(tr).degBinOccupSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).degBinOccup, 2)';
                rawDataByTrialNeural{znType,znNum}(tr).degBinPosSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).degBinPos, 2)';
                rawDataByTrialNeural{znType,znNum}(tr).degBinSpeedSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).degBinSpeed, 2)';
                rawDataByTrialNeural{znType,znNum}(tr).degBinLickRateSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).degBinLickRate, 2)';

                %add spiking data and rate map
                rawDataByTrialNeural{znType,znNum}(tr).degBinSpikeCount = zeros(length(degBinEdges)-1,length(rawDataBySessionNeural.apData));
                rawDataByTrialNeural{znType,znNum}(tr).degBinSpikeCountSmooth = zeros(length(degBinEdges)-1,length(rawDataBySessionNeural.apData));
                rawDataByTrialNeural{znType,znNum}(tr).degBinRateMap = zeros(length(degBinEdges)-1,length(rawDataBySessionNeural.apData));
                for clu = 1:length(rawDataBySessionNeural.apData)
                    spikePosInd = lookup2(rawDataByTrialNeural{znType,znNum}(tr).apData(clu).spikeInds, rawDataByTrialNeural{znType,znNum}(tr).apTime);
                    spikePos = rawDataByTrialNeural{znType,znNum}(tr).currentDeg(spikePosInd);
                    if any(abs(diff(spikePos)) > 300) %for trials that wrap 360
                       spikePos(spikePos<300) = spikePos(spikePos<300) + 360;
                    end
                    rawDataByTrialNeural{znType,znNum}(tr).degBinSpikeCount(:,clu) = histcounts(spikePos,degBinEdges);
                    rawDataByTrialNeural{znType,znNum}(tr).degBinSpikeCountSmooth(:,clu) = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).degBinSpikeCount(:,clu)', 2)';
                    rawDataByTrialNeural{znType,znNum}(tr).degBinRateMap(:,clu) = rawDataByTrialNeural{znType,znNum}(tr).degBinSpikeCountSmooth(:,clu) ./ rawDataByTrialNeural{znType,znNum}(tr).degBinOccupSmooth;
                end%clu

                % % Below code for binning by time works, but currnetly not including time bins because 
                % %  rawDataByTrialNeural is extracted by degrees, leading to many time bins being NaN/0.
                % %  Also, not many spikes, so the several time bins leads to noisy results.
                %TIME%
                % % entryTime = rawDataByTrialNeural{znType,znNum}(tr).apTime(find(rawDataByTrial{znType,znNum}(tr).currentDeg >= floor(rawDataByTrial{znType,znNum}(tr).currentDeg(1))+params.gapBefore,1,'first'));
                % % tBinEdges = entryTime-params.binsize_stime*params.samprate:params.samprate*params.binsize_mstime/1000:entryTime+params.binsize_stime*params.samprate;%in ms bins
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinSizeMS = params.binsize_mstime;
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinSizeSamp = params.samprate*params.binsize_mstime/1000;
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinEdges = tBinEdges;
                % % 
                % % %identify which bins each position data point belongs to
                % % tBinIden = discretize(rawDataByTrialNeural{znType,znNum}(tr).apTime,tBinEdges);
                % % 
                % % %compute occupancy, speed, and licks per spatial bin 
                % % for b = 1:length(tBinEdges)-1
                % %     rawDataByTrialNeural{znType,znNum}(tr).tBinOccup(b) = sum(diff(rawDataByTrialNeural{znType,znNum}(tr).vrTime(tBinIden == b)));
                % %     rawDataByTrialNeural{znType,znNum}(tr).tBinPos(b) = sum(diff(rawDataByTrialNeural{znType,znNum}(tr).currentDeg(tBinIden == b)));
                % %     rawDataByTrialNeural{znType,znNum}(tr).tBinLicked(b) = sum(rawDataByTrialNeural{znType,znNum}(tr).licked(tBinIden == b));
                % % end%bin
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinSpeed = rawDataByTrialNeural{znType,znNum}(tr).tBinPos ./ rawDataByTrialNeural{znType,znNum}(tr).tBinOccup;
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinLickRate = rawDataByTrialNeural{znType,znNum}(tr).tBinLicked ./ rawDataByTrialNeural{znType,znNum}(tr).tBinOccup;
                % % %smooth data
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinOccupSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).tBinOccup, 2)';
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinPosSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).tBinPos, 2)';
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinSpeedSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).tBinSpeed, 2)';
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinLickRateSmooth = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).tBinLickRate, 2)';
                % % 
                % % %add spiking data and rate map
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinSpikeCount = zeros(length(tBinEdges)-1,length(rawDataBySessionNeural.apData));
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinSpikeCountSmooth = zeros(length(tBinEdges)-1,length(rawDataBySessionNeural.apData));
                % % rawDataByTrialNeural{znType,znNum}(tr).tBinRateMap = zeros(length(tBinEdges)-1,length(rawDataBySessionNeural.apData));
                % % for clu = 1:length(rawDataBySessionNeural.apData)
                % %     spikeTimeInd = lookup2(rawDataByTrialNeural{znType,znNum}(tr).apData(clu).spikeInds, rawDataByTrialNeural{znType,znNum}(tr).apTime);
                % %     spikeTime = rawDataByTrialNeural{znType,znNum}(tr).apTime(spikeTimeInd);
                % %     rawDataByTrialNeural{znType,znNum}(tr).tBinSpikeCount(:,clu) = histcounts(spikeTime,tBinEdges);
                % %     rawDataByTrialNeural{znType,znNum}(tr).tBinSpikeCountSmooth(:,clu) = gaussSmooth(rawDataByTrialNeural{znType,znNum}(tr).tBinSpikeCount(:,clu)', 2)';
                % %     rawDataByTrialNeural{znType,znNum}(tr).tBinRateMap(:,clu) = rawDataByTrialNeural{znType,znNum}(tr).tBinSpikeCountSmooth(:,clu) ./ rawDataByTrialNeural{znType,znNum}(tr).tBinOccupSmooth;
                % % end%clu

               %%%%% add data to cell after last zone, like rawDataByTrial struct %%%%%
               rawDataByTrialNeural{znType,size(rawDataByTrial,2)}(znNum + (tr-1)*(size(rawDataByTrial,2)-1)) = rawDataByTrialNeural{znType,znNum}(tr);
            end%tr
        end%znNum
    end%znType

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% Save data %%%%%
    %%%%%%%%%%%%%%%%%%%%%
    filename = [saveNeuralPath '\' 'rawDataByTrialNeural.mat'];
    save(filename, 'rawDataByTrialNeural', '-v7.3');

end%if isfile([virmenSessDataPath '\rawDataByTrial.mat'])

end%fucntion