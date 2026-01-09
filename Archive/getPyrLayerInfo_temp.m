function sessionPyrLayerInfo = getPyrLayerInfo_temp(subj, sessDate, sessNum, params,neuralRawDataPath, saveNeuralPath,  plotPyrLayer)

%% load session data %%
load([neuralRawDataPath '\kilosort4\clusters_allrec.mat'])
load([saveNeuralPath '\rawDataBySessionNeural.mat'])

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
            for tMob = 1:100%length(randMobTime)
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
            for tImmob = 1:100%length(randImTime)
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
[allCluMaxChanSortVal, allCluMaxChanSortInd]  = sort(allPyrMaxChan);

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














%determine windows with ripple amplitude, spike amplitude, and number of pyramidal cells above 2 std
pyrLayers = []; pyrLayerCA1tmp = []; pyrLayerCA3tmp = []; pyrLayerCA1 = []; pyrLayerCA3 = [];
pyrLayers = find(allWPyrCount>=(std(allWPyrCount)*2) & allWRipAmpMed>=(std(allWRipAmpMed)*2)...
    & allWSpAmpMed>=(nanstd(allWSpAmpMed)*2) & abs(allWThetaAmpMed)>=(std(allWThetaAmpMed)*2));
%use channel 80 as cutoff for CA1 (above)
pyrLayerCA1tmp = pyrLayers(pyrLayers>=80);

%%%%% identify CA1 %%%%%
%label the CA1 channel with the highest ripple amplitude as the CA1 pyramidal layer
[~, tmpCA1Ind] = max(filtdata(pyrLayerCA1tmp) - median(filtdata(pyrLayerCA1tmp),2));
pyrLayerCA1 = pyrLayerCA1tmp(tmpCA1Ind);

%%%%% identify CA3 %%%%%
%label the CA3 channel with the highest spike amplitude as the CA3 pyramidal layer
pyrLayerCA3tmp = pyrLayers(pyrLayers<80 & pyrLayers>=pyrLayerCA1-70); %CA3 pyramidal layer always within 70 channels of CA1 pyramidal layer
[~, tmpCA3Ind] = max(allPyrSpAmpsMn(pyrLayerCA3tmp));
pyrLayerCA3 = pyrLayerCA3tmp(tmpCA3Ind);


% % %%%%% save data %%%%%
% % sessionPyrLayerInfo.pyrLayerCA1 = pyrLayerCA1-10:pyrLayerCA1+9;
% % sessionPyrLayerInfo.pyrLayerCA3 = pyrLayerCA3-10:pyrLayerCA3+9;
% % filename = [saveNeuralPath '\' 'sessionPyrLayerInfo.mat'];
% % save(filename, 'sessionPyrLayerInfo', '-v7.3')




% % % % % %%%%% find ripple amplitudes by channel %%%%%
% % % % % %filt params
% % % % % filt_freq = [150 240];
% % % % % data = rawDataBySessionNeural.lfpData(:,randTime:randTime+chparams.Fs*numsecs-1)';
% % % % % %build filt
% % % % % N = 3;
% % % % % Wn = filt_freq ./ (chparams.Fs/2);
% % % % % [b,a] = butter(N,Wn);
% % % % % %filter data in spike frequency range
% % % % % filtdata = filtfilt(b, a, data);
% % % % % filtdata = -1*filtdata; %flip signal so spikes go up
% % % % % if size(filtdata,2)>1
% % % % %     filtdata=filtdata';
% % % % % end
% % % % % 
% % % % % %%%%% loop through windows of channels %%%%%
% % % % % %find median ripple amplitude, median spike amplitude, and number of putative pyramidal cells for each window of n channels
% % % % % numChans = 10;
% % % % % maxRipAmpPerCh = max(filtdata,[],2);
% % % % % medRipAmpPerCh = median(filtdata,2);
% % % % % [~, pyrLayerCA1] = max(abs((maxRipAmpPerCh - medRipAmpPerCh)));
% % % % % for w = 1:pyrLayerCA1-1
% % % % %     thisW = w:w+(numChans-1);
% % % % %     allWRipAmpMed(w) = median(maxRipAmpPerCh(thisW));
% % % % %     allWPyrCount(w) = length(allPyrCluID(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
% % % % %     allWSpAmpMed(w) = median(allPyrSpAmpsMn(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
% % % % % end
% % % % % 
% % % % % [~, tmpCh] = max(allPyrSpAmpsMn);
% % % % % pyrLayerCA3 = allPyrMaxChan(tmpCh);



% % % % %  PROBABLY DELETE
% % % % % allWRipAmpMed = []; allWPyrCount = []; allWSpAmpMed = [];
% % % % % for w = 1:size(data,2)-numChans+1
% % % % %     thisW = w:w+(numChans-1);
% % % % %     allWRipAmpMed(w) = median(maxRipAmpPerCh(thisW));
% % % % %     allWPyrCount(w) = length(allPyrCluID(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
% % % % %     allWSpAmpMed(w) = median(allPyrSpAmpsMn(allPyrMaxChan>=thisW(1) & allPyrMaxChan<=thisW(end)));
% % % % % end
% % % % % %determine windows with ripple amplitude, spike amplitude, and number of pyramidal cells above 2 std
% % % % % pyrLayers = []; pyrLayerCA1tmp = []; pyrLayerCA3tmp = []; pyrLayerCA1 = []; pyrLayerCA3 = [];
% % % % % pyrLayers = find(allWPyrCount>=(std(allWPyrCount)*2) & allWRipAmpMed>=(std(allWRipAmpMed)*2) & allWSpAmpMed>=(nanstd(allWSpAmpMed)*2));
% % % % % %use channel 80 as cutoff between CA3 (above) and CA1 (below)
% % % % % pyrLayerCA3tmp = pyrLayers(pyrLayers<80);
% % % % % pyrLayerCA1tmp = pyrLayers(pyrLayers>80);
% % % % % 
% % % % % %%%%% identify CA1 %%%%%
% % % % % %label the CA1 channel with the highest ripple amplitude as the CA1 pyramidal layer
% % % % % [~,maxCA1PyrCh] = max(allWRipAmpMed(pyrLayerCA1tmp));
% % % % % pyrLayerCA1 = pyrLayerCA1tmp(maxCA1PyrCh);
% % % % % 
% % % % % %%%%% identify CA3 %%%%%
% % % % % %label the CA3 channel with the highest spike amplitude as the CA3 pyramidal layer
% % % % % pyrLayerCA3CluChs = allCluMaxChanSortVal(allPyrMaxChan>=pyrLayerCA3tmp(1) & allPyrMaxChan<=pyrLayerCA3tmp(end));
% % % % % [~,maxCA3PyrCh] = max(allPyrSpAmpsMn(allPyrMaxChan>=pyrLayerCA3tmp(1) & allPyrMaxChan<=pyrLayerCA3tmp(end)));
% % % % % pyrLayerCA3 = pyrLayerCA3CluChs(maxCA3PyrCh);

%% Plot %%
if plotPyrLayer

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% plot depth plot of power %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %includes power by channel and location of putative cells and putative pyramidal layers
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
            %fprintf('\tPlotting ch %d of %d \n', ch, num_chans_to_plot)
            %figure out which y coordinates
            ycoords = [ch-.5 ch-.5 ch+.5 ch+.5];
            tmp_power_mn = mean(all_ch_power(ch,whichf));%mean power over this range of frequencies
            %figure out the color for this value
            myc = ceil((tmp_power_mn-freq(frng).pclim(1))/(freq(frng).pclim(2)-freq(frng).pclim(1)) * length(mycolors(:,1)));
            %cant be less than min or more than max color
            if(myc<1), myc = 1; end
            if(myc>length(mycolors(:,1))), myc = length(mycolors(:,1)); end
            myc = mycolors(myc, :);%the actual RGB values
            %use patch to draw a rectangle
            patch(xcoords, ycoords, myc);
        end%ch

        %plot clusters on top of channels
        if frng == 4
            hold on
            plot(sum(xcoords(1:2))/2, allPyrMaxChan,'^k','MarkerFaceColor', 'k')
            plot(sum(xcoords(1:2))/2+.2, allIntMaxChan,'*c')
            % plot(xcoords(1:2),[pyrLayers; pyrLayers]','--g','LineWidth',2)
            plot(xcoords(1:2),[pyrLayerCA1 pyrLayerCA1],'--g','LineWidth',2)
            plot(xcoords(1:2),[pyrLayerCA3 pyrLayerCA3],'--g','LineWidth',2)
        end%if frng == 4

        %plot details
        xlim([0 5])
        ylim([0 num_chans_to_plot+1])
        yticks(1:2:num_chans_to_plot)
        xticks([1.5 3.5])
        xticklabels({'Theta','Spiking'})
        set(gca, 'FontSize', 7, 'FontName', 'Arial')
        ylabel('Channels', 'FontWeight', 'bold', 'FontSize', 10, 'FontName', 'Arial')
        title(sprintf('Depth Plot: %s_%s_%s',subj,sessDate,sessNum), 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial','Interpreter','none')

    end%frng

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% plot depth plots of raw + filtered traces %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Note: Can add black lines to investigate phase across channels shift
    %plot raw data
    chans = pyrLayerCA1-20:2:pyrLayerCA1+20;%by 4 = 40 microns, so -20 to +20 channels = 600 microns
    %chans = 145-20:2:145+20;
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
    plot_spacer = 20;
    figure; hold on; for ch = 1:length(chans); plot(filtdata(chans(ch),:)'+ch/plot_spacer);end
    title('CA1 Filtered Traces of Random Immobile Time')
    xlabel('Time (s)')
    xticks([1:chparams.Fs-1:chparams.Fs*numsecs-1])
    xticklabels([0 1 2 3])
    ylabel('Channel')
    yticks([1:size(data,2)]/plot_spacer)
    yticklabels([chans])

    % % %%%%% plot power spectra for all channels; dark blue = lower channels, dark red = higher channels
    % % figure
    % % hold on
    % % colors = cbrewer('div','RdBu', size(all_ch_power,1)+3);
    % % colors = flip(colors);
    % % img = arrayfun( @(x) plot(all_ch_power(x,:), 'Color', abs(colors(x+3,:)), 'DisplayName', num2str(x)), 1:size(all_ch_power,1));
    % % xticks(1:20:length(f)); xticklabels(round(f(1:20:end))); xlabel('Frequency (Hz)')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% plot bar chart of spike amplitudes + counts by channel %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure; hold on
    bar(allPyrSpAmpsMn(allCluMaxChanSortInd))%in accending channel order
    plot(allPyrNumSp(allCluMaxChanSortInd)/max(allPyrNumSp)*100,'*k')
    plot(1:length(allPyrSpAmpsMn),ones(1,length(allPyrSpAmpsMn))*median(allPyrSpAmpsMn),'--k')%median number of spikes across clusters
    xticks(1:length(allCluMaxChanSortInd))
    xticklabels(allCluMaxChanSortVal)
    xlabel('Channel')
    ylabel('Mean Spike Amplitude')
    title(sprintf('Pyramidal Cell Info: %s_%s_%s',subj,sessDate,sessNum), 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial','Interpreter','none')
end%if plotPyrLayer

end%function