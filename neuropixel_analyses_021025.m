%Using SGLX_readMeta.m function from https://github.com/jenniferColonell/SpikeGLX_Datafile_Tools/tree/main/MATLAB
experimental = 1;
if experimental == 1
    expStr = 'Neuropixels';
    recFold = 'VR_AnnularReversalTask_JK';
elseif experimental == 0
    expStr = 'Neuropixels_practice';
    recFold = 'Neuropixels_Practice';
end
subj = 'JK15_250530';
sess = '3';

if isfolder('\\ad.gatech.edu')
    pathst = '\\ad.gatech.edu\\bme\\labs\\';
elseif isfolder('Y:')
    pathst = 'Y:\';
end

%% Plot Probe Depth Charts %%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% read lfp data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%read lfp meta
lfp_binName = sprintf('%s_%s_%s_g0_t0.imec0.lf.bin', expStr, subj, sess);
lfp_path = [pathst sprintf('singer\\RawData\\%s\\%s\\%s_%s_%s_g0\\%s_%s_%s_g0_imec0',recFold,subj,expStr,subj,sess,expStr,subj,sess)];

lfp_meta = SGLX_readMeta.ReadMeta(lfp_binName, lfp_path);

%read lfp data
samp0 = 0;
lfp_samprate = str2double(lfp_meta.imSampRate);%Hz
nSamp = str2double(lfp_meta.fileTimeSecs)*lfp_samprate;
lfp_data = SGLX_readMeta.ReadBin(samp0, nSamp, lfp_meta, lfp_binName, lfp_path);

%correct lfp data for gain
lfp_data = SGLX_readMeta.GainCorrectOBX(lfp_data, [1:str2double(lfp_meta.nSavedChans)], lfp_meta);

%%%%%%%%%%%%%%%%%%%%%%
%%%%% parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%
%set chronux params
params.fpass = [1 400];
params.tapers = [2 3];
params.Fs = lfp_samprate;
params.pad = 0;
params.trialave = 0;

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
all_ch_lfp = [];
all_ch_power = [];
mych = 1:str2double(lfp_meta.nSavedChans);%all channels
%mych = 1:round(str2double(lfp_meta.nSavedChans)/2);%bottom half of probe (CA3)
%mych = round(str2double(lfp_meta.nSavedChans)/2):str2double(lfp_meta.nSavedChans)-1;%top half of probe (CA1)
chctr = 0;
for ch = mych
    chctr = chctr + 1;
    %get data for this channel; using random chunk of data
    randtime = 900000;
    numsecs = 3;
    tmp_lfp = lfp_data(ch,randtime:randtime+lfp_samprate*numsecs-1);

    %calculate spectra using Chronux
    %S = frequency bins x windows
    [S, f] = mtspectrumc(tmp_lfp',params);

    %find median power across windows; tmp_power = frequency bins x 1
    df = 2*(size(tmp_lfp,1))*params.tapers(2);%degrees of freedom = 2 * #trials * #tapers
    tmp_power = median(S,2);%use median across events to avoid impact of outlier events

    %save data to all channels struct
    all_ch_lfp(chctr,:) = tmp_lfp;
    all_ch_power(chctr,:) = 10*(log10(tmp_power)-psi(df/2)+log(df/2)); %log transformed and bias corrected and X10 bel-->decibel
end%ch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% plot power for all channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_chans_to_plot = length(mych);
figure
for frng = [1 4]%theta and spiking
    %figure out which frequency bins to use for this range
    whichf = f>=freq(frng).range(1) & f<=freq(frng).range(2);
    %figure out color map
    numcolors = 20;
    mycolors = jet(numcolors);%one color from blue to red for each power unit
    %figure out which x coordinates
    xcoords = [freq(frng).xlim fliplr(freq(frng).xlim)];
    %figure out min and max power across all channels for this spiking range
    maxval = max(all_ch_power(:,whichf),[],'all');
    minval = min(all_ch_power(:,whichf),[],'all');
    freq(frng).pclim = [minval maxval];

    for ch = 1:num_chans_to_plot
        fprintf('\tPlotting ch %d of %d \n', ch, num_chans_to_plot)
        %figure out which y coordinates
        ycoords = [ch-.5 ch-.5 ch+.5 ch+.5];
        tmp_power_mn = mean(all_ch_power(ch,whichf));%mean power over this range of frequencies
        %save power in the spiking range for each channel to determine putative pyramidal layer channel
        if frng == 4 
            all_ch_spike_mn(ch) = tmp_power_mn;
        end

        %figure out the color for this value
        myc = ceil((tmp_power_mn-freq(frng).pclim(1))/(freq(frng).pclim(2)-freq(frng).pclim(1)) * length(mycolors(:,1)));
        %cant be less than min or more than max color
        if(myc<1), myc = 1; end
        if(myc>length(mycolors(:,1))), myc = length(mycolors(:,1)); end
        myc = mycolors(myc, :);%the actual RGB values
        %use patch to draw a rectangle
        patch(xcoords, ycoords, myc);
    end%ch

    %plot details
    xlim([0 5])
    ylim([0 num_chans_to_plot+1])
    yticks(1:2:num_chans_to_plot)
    xticks([1.5 3.5])
    xticklabels({'Theta','Spiking'})
    set(gca, 'FontSize', 7, 'FontName', 'Arial')
    ylabel('Channels', 'FontWeight', 'bold', 'FontSize', 10, 'FontName', 'Arial')
    title(sprintf('Power'), 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial')

    %find putative pyramidal layer channel after removing noisy channels
    if frng == 4
        %remove noisy channels
        noisych_diff = find(abs(diff(all_ch_spike_mn))>std(all_ch_spike_mn)*2);
        if sum(noisych_diff)>0
            for n = 1:length(noisych_diff)
                %find which channel to remove of the 2 channels whose
                % difference is greater than 3 std; done by finding the
                % channel with the largest value after subtracting the
                % median of the three surrounding channels
                tmpnoisych = abs(all_ch_spike_mn(noisych_diff(n)-1:noisych_diff(n)+1) - median(all_ch_spike_mn(noisych_diff(n)-1:noisych_diff(n)+1)));
                [~, tmpnoisych_ind] = max(tmpnoisych);
                if tmpnoisych_ind == 1
                    all_ch_spike_mn(noisych_diff(n)-1) = nan;
                elseif tmpnoisych_ind == 2
                    all_ch_spike_mn(noisych_diff(n)) = nan;
                elseif tmpnoisych_ind == 3
                    all_ch_spike_mn(noisych_diff(n)+1) = nan;
                end%if tmpnoisych_ind == 1
            end%n
        end%if sum(noisych_diff)>0
        %find pyramidal layer now that any noisy channels were removed
        [~,pyr_layer] = max(all_ch_spike_mn,[],'all');
    end%if frng == 4

end%frng

% % %%%%% plot depth plot of raw traces; can add black lines to investigate phase across channels shift 
% % %chans = 1:6:100;
% % chans = pyr_layer-20:4:pyr_layer+20;%by 4 = 40 microns, so -20 to +20 channels = 600 microns
% % figure; hold on; for ch = chans; plot(lfp_data(ch,randtime:randtime+lfp_samprate*numsecs-1)+ch/100'); end
% % xticks([1:lfp_samprate-1:lfp_samprate*numsecs-1])
% % xticklabels([0 1 2 3])
% % xlabel('Time (s)')
% % 
% % %%%%% plot depth plot of spike filtered traces
% % %filt params
% % filt_freq = [150 240];
% % data = lfp_data(chans,randtime:randtime+lfp_samprate*numsecs-1)';
% % samprate = lfp_samprate;
% % %build filt
% % N = 3;
% % Wn = filt_freq ./ (samprate/2);
% % [b,a] = butter(N,Wn);
% % %filter data in spike frequency range
% % filtdata = filtfilt(b, a, data);
% % filtdata = -1*filtdata; %flip signal so spikes go up
% % if size(filtdata,2)>1
% %     filtdata=filtdata';
% % end
% % %plot filtered data
% % plot_spacer = 20;
% % figure; hold on; for lp = 1:size(filtdata,1); plot(filtdata(lp,:)'+lp/plot_spacer);end
% % title('Filtered Traces of Random Time')
% % xlabel('Time (s)')
% % xticks([1:lfp_samprate-1:lfp_samprate*numsecs-1])
% % xticklabels([0 1 2 3])
% % ylabel('Channel')
% % yticks([1:size(data,2)]/plot_spacer)
% % yticklabels([chans])
% % 
% % %%%%% plot power spectra for all channels; dark blue = lower channels, dark red = higher channels
% % figure
% % hold on
% % colors = cbrewer('div','RdBu', size(all_ch_power,1)+3);
% % colors = flip(colors);
% % img = arrayfun( @(x) plot(all_ch_power(x,:), 'Color', abs(colors(x+3,:)), 'DisplayName', num2str(x)), 1:size(all_ch_power,1));
% % xticks(1:20:length(f)); xticklabels(round(f(1:20:end))); xlabel('Frequency (Hz)')


%% Create Virmen Data Struct with LFP Indices for Time %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% read nidaq data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read nidaq meta
nidaq_binName = sprintf('%s_%s_%s_g0_t0.nidq.bin',expStr,subj,sess);
nidaq_path = [pathst sprintf('singer\\RawData\\%s\\%s\\%s_%s_%s_g0',recFold,subj,expStr,subj,sess)];

nidaq_meta = SGLX_readMeta.ReadMeta(nidaq_binName, nidaq_path);

%read nidaq data
samp0 = 0;
nidaq_samprate = str2double(nidaq_meta.niSampRate);%Hz
nSamp = str2double(nidaq_meta.fileTimeSecs)*nidaq_samprate;
nidaq_data = SGLX_readMeta.ReadBin(samp0, nSamp, nidaq_meta, nidaq_binName, nidaq_path);

%correct nidaq data for gain
nidaq_data = SGLX_readMeta.GainCorrectNI(nidaq_data, [1:str2double(nidaq_meta.nSavedChans)], nidaq_meta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% align VR + LFP start times %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find VR start and end values based on TTL pulses received from Virmen by SpikeGLX
nidaq_ttl_start_ind = find(diff(nidaq_data(1,:))>4,1,'first')+1;%TTL pulses sent from Virmen to SpikeGLX are 4.9998 mV
nidaq_ttl_start_time = nidaq_ttl_start_ind/nidaq_samprate;%seconds
nidaq_ttl_end_ind = find(diff(nidaq_data(1,:))>4,1,'last')+1;%TTL pulses sent from Virmen to SpikeGLX are 4.9998 mV
nidaq_ttl_end_time = nidaq_ttl_end_ind/nidaq_samprate;%seconds
nidaq_ttl_duration = nidaq_ttl_end_time - nidaq_ttl_start_time;%seconds
%figure; plot(nidaq_data(1,nidaq_ttl_end_ind-20000:nidaq_ttl_end_ind+5000))%plot last few TTL pulses

%find VR start and end values based on TTL pulses sent by Virmen to SpikeGLX
virmen_path = [pathst sprintf('singer\\Josh\\Behavior\\AnnularFAM\\%s_%s',subj,sess)];
%virmen_path = 'Y:\singer\Josh\Behavior\AnnularFAM\X28_250207_1';
virmen_data = load([virmen_path '\' 'dataWithLickometer.mat']);
trigColumn = find(strcmp(virmen_data.saveDatCopy.dataHeaders,'AutoTrigger'));%11th column is for triggers; first 1 = first trigger
virmen_ttl_start_ind = find(virmen_data.saveDatCopy.data(:,trigColumn)==1,1,'first');
virmen_ttl_start_time = sum([0; diff(virmen_data.saveDatCopy.data(1:virmen_ttl_start_ind,1)) .* 24 * 60 * 60]);%seconds
virmen_ttl_end_ind = find(virmen_data.saveDatCopy.data(:,trigColumn)==1,1,'last');
virmen_ttl_end_time = sum([0; diff(virmen_data.saveDatCopy.data(1:virmen_ttl_end_ind,1)) .* 24 * 60 * 60]);%seconds
virmen_ttl_duration = virmen_ttl_end_time - virmen_ttl_start_time;%seconds
%figure; plot(virmen_data.saveDatCopy.data(virmen_ttl_end_ind-100:virmen_ttl_end_ind+50,11))%plot last few TTL pulses

%ensure the duration of TTL pulses is similar (within a second) for virmen and NIDAQ
if abs(virmen_ttl_duration - nidaq_ttl_duration) > 1%seconds
    sprintf('Check duration of TTL pulses for Virmen vs. NIDAQ.')
end%duration check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% make virmen data struct with lfp indices for time %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lfp_ttl_start_ind = nidaq_ttl_start_time*lfp_samprate;
lfp_ttl_end_ind = nidaq_ttl_end_time*lfp_samprate;
virmen_rec_data = virmen_data.saveDatCopy.data(virmen_ttl_start_ind:virmen_ttl_end_ind,:);
virmen_rec_time_diffs = [0; diff(virmen_rec_data(:,1)) .* 24 * 60 * 60] .* lfp_samprate;
virmen_rec_time_sums = cumsum(virmen_rec_time_diffs);
virmen_rec_data(:,1) = lfp_ttl_start_ind + virmen_rec_time_sums;

%ensure last virmen data index is around last lfp ttl index
if abs( (lfp_ttl_end_ind/lfp_samprate) - (virmen_rec_data(end,1)/lfp_samprate) ) > 1%seconds
    sprintf('Check conversion of virmen times to lfp times.')
end%time check

%% Turn Virmen Data Struct into Lap Struct %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find last index of each lap (include first recorded index) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lap_end_inds = [find(abs(diff(virmen_rec_data(:,2))) >= 300)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% make virmen data into lap struct %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lap = struct([]);
for lp = 1:length(lap_end_inds)+1
    if lp == 1
        lap(lp).data = virmen_rec_data(1:lap_end_inds(lp),:);
    elseif lp > 1 && lp ~= length(lap_end_inds)+1
        lap(lp).data = virmen_rec_data(lap_end_inds(lp-1)+1:lap_end_inds(lp),:);
    elseif lp > 1 && lp == length(lap_end_inds)+1
        lap(lp).data = virmen_rec_data(lap_end_inds(lp-1)+1:end,:);
    end%if lp == 1
end%lp

%% Add RZ and CZ information to Lap Struct %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% use mouse name to determine the original zones %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpstr = strfind(virmen_data.saveDatCopy.sessioninfo, '_');
mouseName = virmen_data.saveDatCopy.sessioninfo(1:tmpstr-1);

%define RZs and CZs
track_RZs = virmen_data.saveDatCopy.thetaReward;
if strcmp(virmen_data.saveDatCopy.trackname, 'TrackA')%original
    track_CZs = virmen_data.saveDatCopy.thetaReward+30;
else%update
    if sum(strcmp(mouseName, {'X28', 'JK5', 'JK6', 'JK7', 'JK8', 'JK9', 'JK10'}))>1
        track_CZs = [50 140 240 320];%original zones
    elseif sum(strcmp(mouseName, {'JK11','JK12'}))>1
        track_CZs = [30 140 240 320];%original zones
    else
        track_CZs = [30 140 240];%original zones
    end%~isempty
end%strcmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% add RZ and CZ virmen data and lfps %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ch_to_use = pyr_layer;
for lp = 1:length(lap)
    %first, RZ
    %include the zone before the RZ and the RZ
    lap_RZ_inds = lap(lp).data(:,2)>track_RZs-10 & lap(lp).data(:,2)<track_RZs+10;%zone = 10 degrees
    lap(lp).RZ = struct([]);
    for z = 1:length(track_RZs)
        %add virmen data to lap struct
        lap(lp).RZ(z).data = lap(lp).data(lap_RZ_inds(:,z),:);
        %note: zone will be empty if not reached

        %add lfp data to lap struct; for now, just putative pyramidal layer
        for ch = 1:length(ch_to_use)
            if ~isempty(lap(lp).RZ(z).data)%reward zone
                lap(lp).RZ(z).lfp(ch,:) = lfp_data(ch_to_use(ch),lap(lp).RZ(z).data(1,1):lap(lp).RZ(z).data(end,1));
                lap(lp).RZ(z).lfp_times(ch,:) = lap(lp).RZ(z).data(1,1):lap(lp).RZ(z).data(end,1);
            end%~isempty
        end%ch
    end%zone

    %next, CZ
    %include the zone before the CZ and the CZ
    lap_CZ_inds = lap(lp).data(:,2)>track_CZs-10 & lap(lp).data(:,2)<track_CZs+10;%zone = 10 degrees
    lap(lp).CZ = struct([]);
    for z = 1:length(track_CZs)
        %add virmen data to lap struct
        lap(lp).CZ(z).data = lap(lp).data(lap_CZ_inds(:,z),:);
        %note: zone will be empty if not reached

        %add lfp data to lap struct; for now, just putative pyramidal layer
        for ch = 1:length(ch_to_use)
            if ~isempty(lap(lp).CZ(z).data)%control zone
                lap(lp).CZ(z).lfp(ch,:) = lfp_data(ch_to_use(ch),lap(lp).CZ(z).data(1,1):lap(lp).CZ(z).data(end,1));
                lap(lp).CZ(z).lfp_times(ch,:) = lap(lp).CZ(z).data(1,1):lap(lp).CZ(z).data(end,1);
            end%~isempty
        end%ch
    end%zone
end%lap

%% RZ and CZ LFP analyses %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% create trial x samples matrix of lfps %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
secbefore = 1;
RZ_lfp_ctr = 0;
CZ_lfp_ctr = 0;

%filters
RZ_to_use = length(track_RZs);
CZ_to_use = length(track_CZs);
corronly = 0;

numsamps = secbefore*lfp_samprate;
all_RZ_lfp = nan((length(lap)*length(track_RZs)), numsamps*2);
all_CZ_lfp = nan((length(lap)*length(track_CZs)), numsamps*2);
for lp = 1:length(lap)
    %first, RZ
    for z = 1:RZ_to_use
        if ~isempty(lap(lp).RZ(z).data)
            %add one to counter
            RZ_lfp_ctr = RZ_lfp_ctr + 1;
            %find indices
            thetaColumn = find(strcmp(virmen_data.saveDatCopy.dataHeaders,'theta'));%2nd column is for theta (degrees)
            RZ_start_ind = find(lap(lp).RZ(z).data(:,thetaColumn)>=track_RZs(z),1,'first');
            RZ_start_time = lap(lp).RZ(z).data(RZ_start_ind,1);
            RZ_lfp_start_ind = find(lap(lp).RZ(z).lfp_times >= RZ_start_time,1,'first');
            
            %put all RZs together in their own matrix
            all_RZ_lfp(RZ_lfp_ctr,:) = lap(lp).RZ(z).lfp(1,(RZ_lfp_start_ind-numsamps):(RZ_lfp_start_ind+numsamps-1));

            %remove if correct only
            numRewardsColumn = find(strcmp(virmen_data.saveDatCopy.dataHeaders,'numRewards'));%5th column is for number of rewards
            if corronly && (lap(lp).RZ(z).data(end,numRewardsColumn) - lap(lp).RZ(z).data(1,numRewardsColumn) == 0)%no rewards received in RZ
                all_RZ_lfp(RZ_lfp_ctr,:) = [];%remove data from matrix
            end%corronly
        end%~isempty
    end%zone

    %next, CZ
    for z = 1:CZ_to_use
        if ~isempty(lap(lp).CZ(z).data)
            %add one to counter
            CZ_lfp_ctr = CZ_lfp_ctr + 1;
            %find indices
            thetaColumn = find(strcmp(virmen_data.saveDatCopy.dataHeaders,'theta'));%2nd column is for theta (degrees)
            CZ_start_ind = find(lap(lp).CZ(z).data(:,thetaColumn)>=track_CZs(z),1,'first');
            CZ_start_time = lap(lp).CZ(z).data(CZ_start_ind,1);
            CZ_lfp_start_ind = find(lap(lp).CZ(z).lfp_times >= CZ_start_time,1,'first');
            
            %put all CZs together in their own matrix
            all_CZ_lfp(CZ_lfp_ctr,:) = lap(lp).CZ(z).lfp(1,(CZ_lfp_start_ind-numsamps):(CZ_lfp_start_ind+numsamps-1));
        end%~isempty
    end%zone
end%lap

%remove rows of NaNs
all_RZ_lfp(isnan(all_RZ_lfp(:,1)),:) = [];
all_CZ_lfp(isnan(all_CZ_lfp(:,1)),:) = [];

% % %%%%% plot raw + spike-filtered data for RZs %%%%%
% % %filt params
% % freq = [150 240];
% % data = all_RZ_lfp';
% % samprate = lfp_samprate;
% % 
% % %build filt
% % N = 3;
% % Wn = freq ./ (samprate/2);
% % [b,a] = butter(N,Wn);
% % 
% % %filter data in SWRs frequency range
% % filtdata = filtfilt(b, a, data);
% % filtdata = -1*filtdata; %flip signal so spikes go up
% % if size(filtdata,2)>1
% %     filtdata=filtdata';
% % end
% % 
% % %plot raw and filtered data
% % figure; hold on; for lp = 1:size(filtdata,1); plot(filtdata(lp,:)'+lp/2.5);end
% % hold on; for lp = 1:size(all_RZ_lfp,1); plot(all_RZ_lfp(lp,:)'+lp/2.5);end
% % title('Raw and Filtered Traces of RZs')
% % xlabel('Samples (samprate = 2500Hz)')
% % ylabel('Trials')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% calculate moving window power spectrogram %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.tapers = [2 3];
params.fpass = [1 250];
movingwin = [0.5 0.1];
[S_RZ,t,f]=mtspecgramc(all_RZ_lfp',movingwin,params);%reward zone
[S_CZ,t,f]=mtspecgramc(all_CZ_lfp',movingwin,params);%control zone

%find median power across events; RZ_power = time x frequency bins
%reward zones
df_RZ = 2*(size(all_RZ_lfp,1))*params.tapers(2);%degrees of freedom = 2 * #trials * #tapers
med_S_RZ = median(S_RZ,3);%use median across events to avoid impact of outlier events
RZ_power = 10*(log10(med_S_RZ)-psi(df_RZ/2)+log(df_RZ/2)); %log transformed and bias corrected and X10 bel-->decibel
%control zones
df_CZ = 2*(size(all_CZ_lfp,1))*params.tapers(2);%degrees of freedom = 2 * #trials * #tapers
med_S_CZ = median(S_CZ,3);%use median across events to avoid impact of outlier events
CZ_power = 10*(log10(med_S_CZ)-psi(df_CZ/2)+log(df_CZ/2)); %log transformed and bias corrected and X10 bel-->decibel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% plot power spectrograms %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reward zones
figure; pcolor(RZ_power'); shading flat; colormap(jet)
yticks(1:5:length(f)); yticklabels(round(f(1:5:end))); ylabel('Frequency (Hz)')
xticks([1 size(RZ_power,1)/2+movingwin(1) size(RZ_power,1)]); xticklabels([round(t(1)-secbefore,2) 0 round(t(end)-secbefore,2)]); xlabel('Time (s)')
title(sprintf('RZ Power using Ch%d', ch_to_use))
%control zones
figure; pcolor(CZ_power'); shading flat; colormap(jet)
yticks(1:5:length(f)); yticklabels(round(f(1:5:end))); ylabel('Frequency (Hz)')
xticks([1 size(CZ_power,1)/2+movingwin(1) size(CZ_power,1)]); xticklabels([round(t(1)-secbefore,2) 0 round(t(end)-secbefore,2)]); xlabel('Time (s)')
title(sprintf('CZ Power using Ch%d', ch_to_use))

% % %plot moving window spectrogram for a channel across a full session
% % thisch = pyr_layer;
% % [S_sess,t,f] = mtspecgramc(lfp_data(thisch,:),movingwin,params);
% % df_sess = 2*(1)*params.tapers(2);%degrees of freedom = 2 * #trials * #tapers
% % sess_power = 10*(log10(S_sess)-psi(df_sess/2)+log(df_sess/2)); %log transformed and bias corrected and X10 bel-->decibel
% % figure; pcolor(sess_power'); shading flat; colormap(jet)
% % yticks(1:5:length(f)); yticklabels(round(f(1:5:end))); ylabel('Frequency (Hz)')
% % title(sprintf('Session Power using Ch%d', ch_to_use))


%% Kilosort4 Pipeline for Neuropixels %%
%NOTE THAT THE KILOSORT SPIKE SORTING PIPELINE COMBINES ALL RECORDINGS
% TOGETHER TO DETECT CLUSTERS THROUGHOUT THE RECORDING DAY, THEN CREATES STRUCTS THAT
% SPLIT THE DATA INTO INDIVIDUAL RECORDING SESSIONS AND THAT COMBINED ALL
% DATA
%currently for single day

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% create single AP .bin file for all recordings on a given day %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steps:
% 1) Copy the .ap.bin file for each recording into a single location
% 2) Open command window and cd to that location
% 3) Use "COPY /B file1.bin + file2.bin file3.bin" 
%   (e.g., COPY /B Neuropixels_practice_JK13_250425_2_g0_t0.imec0.ap.bin + Neuropixels_practice_JK13_250425_3_g0_t0.imec0.ap.bin +...
%   Neuropixels_practice_JK13_250425_4_g0_t0.imec0.ap.bin Neuropixels_practice_JK13_250425_merged_g0_t0.imec0.ap.bin)
% https://stackoverflow.com/questions/53279744/how-to-join-two-binary-files-on-windows

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% feed single AP .bin file into Kilosort4, then manually spike sort the kilosort output using phy2 %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%after feeding the 'merged' ap.bin file into kilosort4, move the file into the kilosort4 folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% collect property information for each recording %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpdirs = dir([pathst sprintf('singer\\RawData\\%s\\%s\\*_%s*',recFold,subj,subj)]);
props = [];
for r = 1:size(tmpdirs,1)
    %read ap meta data
    ap_binName = sprintf([tmpdirs(r).name '_t0.imec0.ap.bin']);
    ap_path = [tmpdirs(r).folder '\' tmpdirs(r).name '\' sprintf([tmpdirs(r).name '_imec0'])];

    ap_meta = SGLX_readMeta.ReadMeta(ap_binName, ap_path);

    %read ap data
    samp0 = 0;
    ap_samprate = str2double(ap_meta.imSampRate);%Hz
    nSamp = (str2double(ap_meta.fileTimeSecs)*ap_samprate);

    %save info
    props.fileNames = tmpdirs;
    props.recLength(r) = nSamp;
    props.sampRate = ap_samprate;
    props.numChan = str2double(ap_meta.nSavedChans);
    props.fileNums(r) = str2double(tmpdirs(r).name(end-3));
end

kilosort_path = [props.fileNames(1).folder '\' 'kilosort4'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% create raw cluster structures using phy2 output %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adapted from makeClusterStructure.m

%%%%% read clustered information %%%%%
spikeInds = readNPY([kilosort_path '\' 'spike_times.npy']);
spikeID = readNPY([kilosort_path '\' 'spike_clusters.npy']);
[clusterID, clusterGroup] = readClusterGroupsCSV([kilosort_path '\' 'cluster_group.tsv']);
templates = readNPY([kilosort_path '\' 'templates.npy']);
spikeTemplates = readNPY([kilosort_path '\' 'spike_templates.npy']);
channelMap = readNPY([kilosort_path '\' 'channel_map.npy']);
channelPositions = readNPY([kilosort_path '\' 'channel_positions.npy']);
clusterParams = loadParamsPy([kilosort_path '\' 'params.py']);    

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
        else
            tempSpikeInds{clu} = tempSpikeInds{clu} - props.recLength(r-1); %align to start time of this recording
            tmprawclusters(clu).spikeInds = [];
        end
        
        tmprawclusters(clu).spikeInds = tempSpikeInds{clu}(tempSpikeInds{clu} <= props.recLength(r));
        
        if r <  size(props.fileNames,1)
            tempSpikeInds{clu} = tempSpikeInds{clu}(tempSpikeInds{clu} > props.recLength(r));
        end
        tmprawclusters(clu).file = props.fileNums(r);
    end
    %output = rawclusters{session}(cluster)
    rawclusters{r} = tmprawclusters;
end

%make structure with all spike times from all recordings
rawclusters_allrec = struct('ID', num2cell(goodUnits), ...
    'spikeInds', repmat({[]}, 1, length(goodUnits)),...
    'sampRate', num2cell(props.sampRate*ones(1, length(goodUnits))), ...
    'maxChan', num2cell(unitMaxChan'), 'info', repmat({'all files. pre quality control metrics'}, 1, length(goodUnits)), ...
    'numShanks',num2cell(numShanks*ones(1, length(goodUnits))), 'brainReg', repmat({brainReg}, 1, length(goodUnits)), ...
    'numChan', num2cell(props.numChan*ones(1,length(goodUnits))));
[rawclusters_allrec(1:length(goodUnits)).index] = deal(subj); 
[rawclusters_allrec(1:length(goodUnits)).files] = deal(props.fileNums); 

for clu = 1:length(goodUnits)
    tempSpikeInds{clu} = spikeInds(spikeID == goodUnits(clu));
    tempSpikeInds{clu} = double(tempSpikeInds{clu});
    rawclusters_allrec(clu).spikeInds = tempSpikeInds{clu};
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
mmf = memmapfile([kilosort_path '\' dir([kilosort_path '\' sprintf('%s_%s*',expStr,subj)]).name], 'Format', {dataType, [nChInFile nSamp], 'x'});

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
    metrics(clu).index = subj;
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
    clusters(r).index = [subj '_' num2str(fileNum)];
    numgood = 1;

    for clu = 1:length(good_final)
        if good_final(clu)
            clusters(r).data(numgood).ID = rawclusters{r}(clu).ID;
            clusters(r).data(numgood).maxChan = rawclusters{r}(clu).maxChan;
            clusters(r).data(numgood).spikeInds = rawclusters{r}(clu).spikeInds;
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

%save the good clusters + metrics
save([pathst sprintf('singer\\RawData\\%s\\%s\\kilosort4\\clusters.mat',recFold,subj)],'clusters','-v7.3')
save([pathst sprintf('singer\\RawData\\%s\\%s\\kilosort4\\clustermetrics.mat',recFold,subj)],'clustermetrics','-v7.3')
save([pathst sprintf('singer\\RawData\\%s\\%s\\kilosort4\\clusters_allrec.mat',recFold,subj)],'clusters_allrec','-v7.3')

%% Cell type classification via CellExplorer %%
%adapted from CellExplorerSingerLabTesting.m
%Note: CellExplorer uses Kilosort output, and thus uses data from all sessions
% on a given day to classify putative pyramidal cells and interneurons.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Create and Validate Session Struct %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adapted from NovelTaskTemplate.m
session = [];
%%%%% add general metadata %%%%%
session.general.basePath =  kilosort_path; % Full path
session.general.name = subj; % Session name / basename
session.general.version = 5; % Metadata version
session.general.sessionType = 'Acute'; % Type of recording: Chronic, Acute, Unknown
session.general.sessionName = session.general.name;

%%%%% add animal data %%%%%
session.animal.name = session.general.name; % Animal name is inferred from the data path
session.animal.sex = 'Female'; % Male, Female, Unknown
session.animal.species = 'Mouse'; % Mouse, Rat
session.animal.strain = 'APOE';
session.animal.geneticLine = '';

%%%%% add extracellular data %%%%%
session.extracellular.fileName = dir([kilosort_path '\' sprintf('%s_%s*',expStr,subj)]).name;
session.extracellular.leastSignificantBit = 0.1950; %for Neuropixels
session.extracellular.probeDepths = 0;
session.extracellular.precision = 'int16';
session.extracellular.sr = ap_samprate;
session.extracellular.nChannels = length(channelMap);
session.extracellular.nElectrodeGroups = 1;
session.extracellular.electrodeGroups.channels = {1:session.extracellular.nChannels};
session.extracellular.nSpikeGroups = 1;
session.extracellular.spikeGroups.channels = {1:session.extracellular.nChannels};
session.extracellular.chanCoords.x = channelPositions(:,1);
session.extracellular.chanCoords.y = channelPositions(:,2);
session.extracellular.chanCoords.source = 'Kilosort';
session.extracellular.chanCoords.verticalSpacing = 20;
session.extracellular.chanCoords.layout = 'staggered';% Probe layout: linear,staggered,poly2,edge,poly3,poly5

%%%%% add spikeSorting data %%%%%
session.spikeSorting{1}.relativePath = '';
session.spikeSorting{1}.format = 'Phy';
session.spikeSorting{1}.method = 'Kilosort';
session.spikeSorting{1}.channels = [];
session.spikeSorting{1}.manuallyCurated = 1;
session.spikeSorting{1}.notes = '';

%add brainRegions data
session.brainRegions.HIP.channels = 1:session.extracellular.nChannels;
session.brainRegions.HIP.electrodeGroups = 1;

%%%%% add epochs data %%%%%
session.epochs{1}.name = session.general.name;
session.epochs{1}.startTime = 0;

%%%%% validate session struct %%%%%
validateSessionStruct(session);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Run the cell metrics pipeline %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exclude_metrics = {'monoSynaptic_connections','spatial_metrics','event_metrics','manipulation_metrics'};
ProcessCellMetrics('session', session, 'excludeMetrics', exclude_metrics,'forceReload', false, 'keepCellClassification', ~false);
%Output notes: This function produces the following stucts: cell_metrics, noiseLevel, session, and spikes. 
% The function includes only clusters manually spike sorted as "good". Clusters, clusters_allrec, and
%clustermetrics take "good" clusters from manual spike sorting "good" and removes clusters that do not meet lab-determined criteria.
%Thus, clusters, clusters_allrec, and clustermetrics include fewer clusters than cell_metrics, noiseLevel, session, and spikes structs, and
%clusters_allrec.ID can be used to index those structs ( e.g., spikes.ts(ismember(spikes.cluID, [clusters_allrec.ID])) ).
%  For spike times, cell_metrics.spikes.times == spikes.times in seconds (so multiplying by 30k will convert to samples),
%and spikes.ts are already in samples.














    
% % Below is an uncompleted, streamlined version of Abby's celltype
% classification pipeline that she used instead of CellExplorer
% % 
% % %% Cell type classification %%
% % %adapted from cf_celltypeclassifier.m and subfunctions
% % %currently for single day
% % 
% % cellTypeProps = [];
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%% get firing rate data %%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % cellTypeProps.dayindex = clustermetrics(1).index;
% % cellTypeProps.meanFR  = [clustermetrics.firingrate]; %ALP 12/14/2022
% % cellTypeProps.waveforms = clustermetrics2mat(clustermetrics, 'WF', 'mn');
% % 
% % % % %%%%% plot %%%%%
% % % % figure
% % % % edges = 0:0.5:30;
% % % % histogram(cellTypeProps.meanFR, edges)
% % % % ylabel('Number of Clusters'); xlabel('Mean FR');
% % % % title(['FR distribution - ', subj], 'Interpreter', 'none');
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%% get spikewidth data %%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %adapted from getspikewidthforcelltypeclassification_K2.m and calcSpikewidth_K2.m
% % for clu = 1:length(clustermetrics)
% %         if ~isnan(clustermetrics(clu).WF.mn(1))
% %             %caluclate peak to trough
% %             resamplefactor = round(16/(props.sampRate/10000)); %8 for 20kHz, 5 for 30kHz (rounding 5.333), should optimize ALP 7/25
% %             waveformresampled = resample(clustermetrics(clu).WF.mn,resamplefactor,1);
% %             % waveformresampled = waveformresampled((resamplefactor*15):(60*resamplefactor)); %what do these numbers mean??
% %             waveformresampled = waveformresampled(round((0.75/1000)*resamplefactor*props.sampRate):round(((3/1000)*resamplefactor*props.sampRate))); %for variable samprate
% %             temp = (waveformresampled - mean(waveformresampled));
% %             waveformresampled = temp/abs(min(temp));
% % 
% %             %find peak2trough
% %             peakIdx = find(waveformresampled == min(waveformresampled));
% %             offsetIdx = peakIdx + 2; %don't start looking until 3 samples away from peak
% %             negslope = find(diff(waveformresampled(offsetIdx:end))<0);
% %             diffWF = diff(waveformresampled(offsetIdx:end));
% %             if isempty(negslope)
% %                 negslope = find(diffWF == min(diffWF));
% %             end
% %             troughIdx = negslope(1)+offsetIdx;
% %             SW = 1000*((abs(troughIdx-peakIdx)/(resamplefactor*props.sampRate))); %1000 to convert to msec, resamp*samprate to correct for sampling
% %             waveform = waveformresampled;
% % 
% %             %find peak2troughsimple
% %             troughIdx2 = find(waveformresampled == max(waveformresampled(offsetIdx:end)));
% %             SW2 = 1000*(abs(troughIdx2-peakIdx)/(resamplefactor*props.sampRate));
% % 
% %             %store the info
% %             cellTypeProps.spikewidth.peak2troughDiff(clu) = SW;
% %             cellTypeProps.spikewidth.peakIdxDiff(clu) = peakIdx;
% %             cellTypeProps.spikewidth.troughIdxDiff(clu) = troughIdx;
% %             cellTypeProps.spikewidth.fullWFDiff(clu,:) = waveform;
% %             cellTypeProps.spikewidth.SW2(clu) = SW2;
% %         else
% %             cellTypeProps.spikewidth.peak2troughDiff(clu) = nan;
% %             cellTypeProps.spikewidth.peakIdxDiff(clu) = nan;
% %             cellTypeProps.spikewidth.troughIdxDiff(clu) = nan;
% %             cellTypeProps.spikewidth.fullWFDiff(clu) = nan;
% %             cellTypeProps.spikewidth.SW2(clu) = nan;
% %         end
% % 
% % end
% % 
% % %remove bad waveforms/clusters
% % prompt = {['Enter MatIdx of bad units for:' 'JK13' ' ' '250417']};
% % prompttitle = 'Exclude bad WFs for exclusion';
% % definput = {''}; dims = [1 35];
% % answer = inputdlg(prompt,prompttitle,dims,definput);
% % badWFs = str2num(answer{1});
% % 
% % cellTypeProps.spikewidth.peak2troughDiff(badWFs) = nan;
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%% get autocorrelogram data %%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %adapted from getautocorrelogramforcelltypeclassification_K2.m and calcAutocorr_K2_celltypes.m
% % 
% % %%%%% find autocorrelogram and their mean for each cluster %%%%%
% % for clu = 1:size(clusters_allrec,2)
% %     %make the spike train - using indices instead of times
% %     stepsize = 5 * props.sampRate / 1000; %number of samples for 5ms binsize
% %     spiketrainedges = 0:stepsize:sum(props.recLength); %5ms bins for length of all recordings this day
% %     spiketrain{clu} = histc(clusters_allrec(clu).spikeInds', spiketrainedges);
% % 
% %     %get the autocorr
% %     lag_num = 50 * props.sampRate / 1000; %number of samples for 50ms
% %     lag = lag_num/stepsize; %in bins
% %     autocorr{clu} = xcorr(spiketrain{clu},lag);
% %     autocorr{clu}(lag+1) = 0;
% % 
% %     %eliminate ones with not enough spikes
% %     if max(autocorr{clu}) < 10
% %         autocorr{clu} = nan(1,length(autocorr{clu}));
% %     end
% % 
% %     %get the first moment of the autocorr
% %     sampN = stepsize:stepsize:lag_num;
% %     centerofmass{clu} = (sum(autocorr{clu}(lag+2:end).*sampN)/sum(autocorr{clu}(lag+2:end)))/stepsize;
% % end
% % 
% % %save to cellTypeProps
% % cellTypeProps.autocorr_mean = centerofmass;
% % cellTypeProps.autocorr = autocorr;
% % 
% % %%%%% plot %%%%%
% % %plot the individual autocorrelograms and their respective WFs
% % %CAUTION! LOTS OF PLOTS!
% % for clu = 1:5%:size(clusters_allrec,2)
% % 
% %     binsize = 5;
% %     time = (binsize*(length(cell2mat(cellTypeProps.autocorr(clu)))-1))/2;
% %     figure; clf;
% %     subplot(1,2,1)
% %     bar(-time:binsize:time, cell2mat(cellTypeProps.autocorr(clu))); hold on;
% %     plot(cell2mat(cellTypeProps.autocorr_mean(clu))*5,2,'*g')
% %     title('Autocorrelogram with center of mass (green)')
% %     xlim([-time time])
% %     hold on;
% %     subplot(1,2,2)
% %     plot(cellTypeProps.waveforms')
% %     title(['Waveform - cluster ' num2str(clu) ' ', subj], 'Interpreter', 'none');
% %     xlim([0 size(cellTypeProps.waveforms,2)])
% % 
% % end
% % 
% % % % %plot distribution of autocorrelogram means
% % % % figure; clf; hold on;
% % % % edges = 0:0.2:10;
% % % % histogram(cell2mat(cellTypeProps.autocorr_mean), edges);
% % % % ylabel('Number of Clusters'); xlabel('Autocorrelogram Mean')
% % % % 
% % % % %plot distribution of autocorrelogram means, spikewidth, and mean firing
% % % % % rate; adapted from plotparametersforcelltypeclassification_K2.m
% % % % figure
% % % % plot3(cell2mat(cellTypeProps.autocorr_mean), cellTypeProps.spikewidth.peak2troughDiff, cellTypeProps.meanFR, 'o')
% % % % xlabel('Mean ac (ms)')
% % % % ylabel('Spikewidth - peak2trough diff')
% % % % zlabel('Mean FR (Hz)')