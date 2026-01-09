function [rawDataBySessionNeural, rawDataByLapNeural, rawDataByTrialNeural] = getNeuralStructs_linearJLK(subj,sessDate,sessNum,virmenSessDataPath, neuralRawDataPath, saveNeuralPath)
%adapted from getRestSessionStats_linearJLK

%% Create Raw Data Struct for Whole Session with LFP and AP Times and Data %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load raw virmen data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawDataBySession = load([virmenSessDataPath '\rawDataBySession.mat']);%whole session
rawDataBySession = rawDataBySession.rawDataBySession;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% read nidaq data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read nidaq meta
nidaq_binName = sprintf('Neuropixels_%s_%s_%s_g0_t0.nidq.bin',subj,sessDate,sessNum);
nidaq_path = [neuralRawDataPath '\' sprintf('Neuropixels_%s_%s_%s_g0',subj,sessDate,sessNum)];
nidaq_meta = SGLX_readMeta.ReadMeta(nidaq_binName, nidaq_path);

%read nidaq data
samp0 = 0;
nidaq_samprate = str2double(nidaq_meta.niSampRate);%Hz
nSamp = str2double(nidaq_meta.fileTimeSecs)*nidaq_samprate;
nidaq_data = SGLX_readMeta.ReadBin(samp0, nSamp, nidaq_meta, nidaq_binName, nidaq_path);

%correct nidaq data for gain
nidaq_data = SGLX_readMeta.GainCorrectNI(nidaq_data, [1:str2double(nidaq_meta.nSavedChans)], nidaq_meta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% align VR + nidaq start times %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find VR start and end values based on TTL pulses received from Virmen by SpikeGLX
nidaq_ttl_start_ind = find(diff(nidaq_data(1,:))>4,1,'first')+1;%TTL pulses sent from Virmen to SpikeGLX are 4.9998 mV
nidaq_ttl_start_time = nidaq_ttl_start_ind/nidaq_samprate;%seconds
nidaq_ttl_end_ind = find(diff(nidaq_data(1,:))>4,1,'last')+1;%TTL pulses sent from Virmen to SpikeGLX are 4.9998 mV
nidaq_ttl_end_time = nidaq_ttl_end_ind/nidaq_samprate;%seconds
nidaq_ttl_duration = nidaq_ttl_end_time - nidaq_ttl_start_time;%seconds
%figure; plot(nidaq_data(1,nidaq_ttl_end_ind-20000:nidaq_ttl_end_ind+5000))%plot last few TTL pulses

%find VR start and end values based on TTL pulses sent by Virmen to SpikeGLX
virmen_ttl_start_ind = find(rawDataBySession.ttl==1,1,'first');
virmen_ttl_start_time = sum([0; diff(rawDataBySession.vrTime(1:virmen_ttl_start_ind))]);%seconds
virmen_ttl_end_ind = find(rawDataBySession.ttl==1,1,'last');
virmen_ttl_end_time = sum([0; diff(rawDataBySession.vrTime(1:virmen_ttl_end_ind))]);%seconds
virmen_ttl_duration = virmen_ttl_end_time - virmen_ttl_start_time;%seconds
%figure; plot(rawDataBySession.ttl(virmen_ttl_end_ind-100:virmen_ttl_end_ind+50))%plot last few TTL pulses

%ensure the duration of TTL pulses is similar (within a second) for virmen and NIDAQ
if abs(virmen_ttl_duration - nidaq_ttl_duration) > 1%seconds
    sprintf('Check duration of TTL pulses for Virmen vs. NIDAQ.')
end%duration check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% determine lfp and ap times for virmen session and extract neural data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Determine virmen recodring times %%%%%
virmen_rec_time = rawDataBySession.vrTime(virmen_ttl_start_ind:virmen_ttl_end_ind,:);
virmen_rec_time_diffs = [0; diff(virmen_rec_time)];
virmen_rec_time_sums = cumsum(virmen_rec_time_diffs);

%%%%% LFPs %%%%%
%read lfp meta
lfp_binName = sprintf('Neuropixels_%s_%s_%s_g0_t0.imec0.lf.bin',subj,sessDate,sessNum);
lfp_path = [neuralRawDataPath '\' sprintf('Neuropixels_%s_%s_%s_g0',subj,sessDate,sessNum) '\' sprintf('Neuropixels_%s_%s_%s_g0_imec0',subj,sessDate,sessNum)];
lfp_meta = SGLX_readMeta.ReadMeta(lfp_binName, lfp_path);

%read lfp data
samp0 = 0;
lfp_samprate = str2double(lfp_meta.imSampRate);%Hz
nSamp = str2double(lfp_meta.fileTimeSecs)*lfp_samprate;
lfp_data = SGLX_readMeta.ReadBin(samp0, nSamp, lfp_meta, lfp_binName, lfp_path);

%correct lfp data for gain
lfp_data = SGLX_readMeta.GainCorrectOBX(lfp_data, [1:str2double(lfp_meta.nSavedChans)], lfp_meta);

%determine lfp start/end times
lfp_ttl_start_ind = ceil(nidaq_ttl_start_time*lfp_samprate);
lfp_ttl_end_ind = floor(nidaq_ttl_end_time*lfp_samprate);
virmen_rec_time_lfp = lfp_ttl_start_ind + floor(virmen_rec_time_sums*lfp_samprate);

%ensure last virmen data index is around last lfp ttl index
if abs( (lfp_ttl_end_ind/lfp_samprate) - (virmen_rec_time_lfp(end)/lfp_samprate) ) > 1%seconds
    sprintf('Check conversion of virmen times to lfp times.')
end%time check

%%%%% APs %%%%%
%load session clusters struct
load([neuralRawDataPath '\kilosort4\clusters.mat'])
load([neuralRawDataPath '\kilosort4\' sprintf('%s_%s.cell_metrics.cellinfo.mat', subj, sessDate)])

%extract ap data
for cl = 1:length(clusters); clustersRow(cl) = strcmp(clusters(cl).index,sprintf('%s_%s_%s',subj,sessDate,sessNum));end%cl
ap_samprate = clusters(clustersRow).samprate;
ap_data = clusters(clustersRow).data;

%determine ap start/end times
ap_ttl_start_ind = ceil(nidaq_ttl_start_time*ap_samprate);
ap_ttl_end_ind = floor(nidaq_ttl_end_time*ap_samprate);
virmen_rec_time_ap = ap_ttl_start_ind + floor(virmen_rec_time_sums*ap_samprate);

%ensure last virmen data index is around last ap ttl index
if abs( (ap_ttl_end_ind/ap_samprate) - (virmen_rec_time_ap(end)/ap_samprate) ) > 1%seconds
    sprintf('Check conversion of virmen times to lfp times.')
end%time check

%remove spikes outside of virmen session
for clu = 1:length(ap_data)
    ap_data(clu).spikeInds = ap_data(clu).spikeInds(ap_data(clu).spikeInds >= ap_ttl_start_ind & ap_data(clu).spikeInds <= ap_ttl_end_ind);
end%clu

%add putative cell type
for clu = 1:length(ap_data)
    ap_data(clu).putativeCellType = cell_metrics.putativeCellType{cell_metrics.cluID == ap_data(clu).ID};    
end%clu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Create sturct and save data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawDataBySessionNeural = [];
rawDataBySessionNeural.vrTime = virmen_rec_time;
rawDataBySessionNeural.lfpTime = virmen_rec_time_lfp;
rawDataBySessionNeural.apTime = virmen_rec_time_ap;
rawDataBySessionNeural.lfpData = lfp_data(:,lfp_ttl_start_ind:lfp_ttl_end_ind);
rawDataBySessionNeural.apData = ap_data;
rawDataBySessionNeural.currentDeg = rawDataBySession.currentDeg(virmen_ttl_start_ind:virmen_ttl_end_ind);
rawDataBySessionNeural.speed = [0; diff(rawDataBySessionNeural.currentDeg) ./ diff(rawDataBySessionNeural.vrTime)];
rawDataBySessionNeural.speedSmooth = gaussSmooth(rawDataBySessionNeural.speed', 5)';
rawDataBySessionNeural.rewarded = [0; diff(rawDataBySession.rewards(virmen_ttl_start_ind:virmen_ttl_end_ind))];
rawDataBySessionNeural.licked = [0; diff(rawDataBySession.licks(virmen_ttl_start_ind:virmen_ttl_end_ind))];
rawDataBySessionNeural.currentZone = rawDataBySession.currentZone(virmen_ttl_start_ind:virmen_ttl_end_ind);

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
    
    %loop through laps
    rawDataByLapNeural = [];
    for ii = 1:length(rawDataByLap)
    
        %index whole session neural struct using virmen lap times
        virmen_lap_start_ind = find(rawDataBySessionNeural.vrTime >= rawDataByLap(ii).vrTime(1),1,'first');
        virmen_lap_end_ind = find(rawDataBySessionNeural.vrTime <= rawDataByLap(ii).vrTime(end),1,'last');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% Create sturct %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        rawDataByLapNeural(ii).vrTime = rawDataBySessionNeural.vrTime(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).lfpTime = rawDataBySessionNeural.lfpTime(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).apTime = rawDataBySessionNeural.apTime(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).lfpData = rawDataBySessionNeural.lfpData(:,rawDataByLapNeural(ii).lfpTime(1):rawDataByLapNeural(ii).lfpTime(end));
        for clu = 1:length(rawDataBySessionNeural.apData)
            rawDataByLapNeural(ii).apData(clu).ID = rawDataBySessionNeural.apData(clu).ID;
            rawDataByLapNeural(ii).apData(clu).maxChan = rawDataBySessionNeural.apData(clu).maxChan;
            rawDataByLapNeural(ii).apData(clu).spikeInds = rawDataBySessionNeural.apData(clu).spikeInds(rawDataBySessionNeural.apData(clu).spikeInds >= rawDataByLapNeural(ii).apTime(1) & rawDataBySessionNeural.apData(clu).spikeInds <= rawDataByLapNeural(ii).apTime(end));
            rawDataByLapNeural(ii).apData(clu).putativeCellType = rawDataBySessionNeural.apData(clu).putativeCellType;
        end%clu
        rawDataByLapNeural(ii).currentDeg = rawDataBySessionNeural.currentDeg(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).speed = rawDataBySessionNeural.speed(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).speedSmooth = rawDataBySessionNeural.speedSmooth(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).rewarded = rawDataBySessionNeural.rewarded(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).licked = rawDataBySessionNeural.licked(virmen_lap_start_ind:virmen_lap_end_ind);
        rawDataByLapNeural(ii).currentZone = rawDataBySessionNeural.currentZone(virmen_lap_start_ind:virmen_lap_end_ind);
    
    end%ii
    
    %%%%%%%%%%%%%%%%%%%%%
    %%%%% Save data %%%%%
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
    rawDataByTrial = load([virmenSessDataPath '\rawDataByTrial.mat']);%laps
    rawDataByTrial = rawDataByTrial.rawDataByTrial;
    
    %loop through trials for each zone
    rawDataByTrialNeural = cell(size(rawDataByTrial));
    for znType = 1:size(rawDataByTrial,1)
        for znNum = 1:size(rawDataByTrial,2)
            for tr = 1:length(rawDataByTrial{znType,znNum})
            
                %index whole session neural struct using virmen trial times
                virmen_trial_start_ind = find(rawDataBySessionNeural.vrTime >= rawDataByTrial{znType,znNum}(tr).vrTime(1),1,'first');
                virmen_trial_end_ind = find(rawDataBySessionNeural.vrTime <= rawDataByTrial{znType,znNum}(tr).vrTime(end),1,'last');

                %%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%% Create sturct %%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%
                rawDataByTrialNeural{znType,znNum}(tr).vrTime = rawDataBySessionNeural.vrTime(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).lfpTime = rawDataBySessionNeural.lfpTime(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).apTime = rawDataBySessionNeural.apTime(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).lfpData = rawDataBySessionNeural.lfpData(:,rawDataByTrialNeural{znType,znNum}(tr).lfpTime(1):rawDataByTrialNeural{znType,znNum}(tr).lfpTime(end));
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
                rawDataByTrialNeural{znType,znNum}(tr).licked = rawDataBySessionNeural.licked(virmen_trial_start_ind:virmen_trial_end_ind);
                rawDataByTrialNeural{znType,znNum}(tr).currentZone = rawDataBySessionNeural.currentZone(virmen_trial_start_ind:virmen_trial_end_ind);

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