function rhd2mat_tempbin_DC(datafolder, destinationfolder, index, params)
%rhd2mat_temp bin Extract into files into standard singer lab format. 
%   rh2mat_tempbin takes each rhd file, extracts the data, and concatenates
%   it using a series of temporary bin files. after the data has been saved
%   as raweeg20k.mat, raweeg.mat, as well as the various TTL and ADC data,
%   the temporary files are deleted. 
%   I recommend running extraction to save to local dir vs. server dir for
%   speed. 

% INPUTS
% datafolder --
% destination folder -- folder where extracted matlab files will be saved,
% index -- [recordingfile#] removed subj+date because unnecessary for only
% creating intaninfo
% params
%
% based on extracteegTTLpulses2 ASinger 8/11/15
% ALP 9/17/2020

%% Map the directory, initialize structures and file locations
%%% initialize destination folder
if ~exist(destinationfolder)
    mkdir(destinationfolder); 
end

%%% find the right rhd files
files=dir(datafolder);
for j = 1:length(index)
isrightfile=zeros(length(files),1);
    indrightfile=[];
for f=1:length(isrightfile)
    filename = strsplit(files(f).name, '_');
    if strcmp(filename(1), strcat('recording', num2str(index(j)))) %get the indices of the correct files
        isrightfile(f)=1;
        indrightfile=[indrightfile; f];
    end
end

numfiles=sum(isrightfile);
disp('Loading Raw Data')

%% initialize temporary files, use w permission to wipe anything inside
%amplifier data
tempdir = [destinationfolder, 'tempfiles\'];
if ~isfolder(tempdir)
    mkdir(tempdir);
end

for ch = 1:length(cell2mat(params.savechnum))
    chfilename = [tempdir, 'tempdata', num2str(index(j)), '_ch', num2str(ch-1), '.bin'];
    ampfid{ch} = fopen(chfilename, 'w+');
end

% other data types
digfilename = [tempdir, 'tempdig', num2str(index(j)), '.bin'];
adcfilename = [tempdir, 'tempadc', num2str(index(j)), '.bin'];
digfid = fopen(digfilename, 'w+');
adcfid = fopen(adcfilename, 'w+');

fclose('all'); 
clear ampfid adcfid digfid

%% open files for appending
%amplifier data
%I don't know if I need to do this if its already open as w+...???
%test this!
for ch = 1:length(cell2mat(params.savechnum))
    chfilename = [tempdir, 'tempdata', num2str(index(j)), '_ch', num2str(ch-1), '.bin'];
    ampfid{ch} = fopen(chfilename, 'a+');
end

% other data types
digfid = fopen(digfilename, 'a+');
adcfid = fopen(adcfilename, 'a+');

%% extract the data
t = tic;

for i=1:sum(isrightfile) %for all files in the recording   
    fprintf(1, strcat('Loading...', num2str(i), ' of...', num2str(numfiles), '\n'));   
    
    %%% extract data from .rhd files
    datafile = strcat(datafolder, '\', files(indrightfile(i)).name);
    %commented out for now (next 6 lines), probably want to run all files
    %through gap check. alp 3/16/2021
    %intandata{i} = read_Intan_RHD2000_fileALP(datafile); %read intan data from rhd file
%     if index(1,2) == 201103
%         intandata{i} = read_Intan_RHD2000_fileALP_interpgaps(datafile);
%     else
%         intandata{i} = read_Intan_RHD2000_fileALP(datafile);
%     end
% - end comment 3/16/21
    
    [intandata{i}, gapInfo] = read_Intan_RHD2000_fileALP_interpgaps(datafile); 
    
    %add some helpful info to the structure to save with the raw data
    %structure down below
    if gapInfo.isGap 
        %what sample did this happen at? 
        iGap = gapInfo.iGapStart;
        iGap = iGap + sum(datapointsrec); %adjust to get the sample in the overall file 
        saveGapInfo.isGap = 1; 
        saveGapInfo.nMissingSamp = gapInfo.nMissingSamp;
        saveGapInfo.file = gapInfo.file; 
        saveGapInfo.iGapFullFile = iGap; 
        saveGapInfo.refSampRate = intandata{i}.fileinfo.amplifier_sample_rate; 
    end
    
    %%% get total data points in the recording and other stuff 
    datapointsrec(i) = size(intandata{i}.amplifier_data(1,:),2);
    nAmpChannels = size(intandata{i}.amplifier_data,1); %rows are channels
    nDigInChannels = size(intandata{i}.board_dig_in_data,1);
    nADCChannels = size(intandata{i}.board_adc_data,1);
    
    %%% check channel number
    if ~isequal(nAmpChannels, length(cell2mat(params.savechnum)))
        error('channel number specified in params does not match channels found in file')
    end
            
    %%% write to the txtfile
    for ch = 1:nAmpChannels
        fwrite(ampfid{ch}, intandata{i}.amplifier_data(ch,:), 'double');
    end
    fwrite(digfid, intandata{i}.board_dig_in_data, 'double');
    fwrite(adcfid, intandata{i}.board_adc_data, 'double');
    
    %%% make last 100 pts vector to check for correct reading
    if i == sum(isrightfile)
        last100 = intandata{i}.amplifier_data(1,end-100:end); %channel 0 of port A
        diglast100 = intandata{i}.board_dig_in_data(1,end-100:end); 
        adclast100 = intandata{i}.board_adc_data(1,end-100:end); 
    end
    
    %%% check to make sure these files actually go together
    if i>1
        if (intandata{i}.t_amplifier(1) - intandata{i-1}.t_amplifier(end))<=0.01
            intandata{i-1} = []; %clear out data
        else
            error('These files do not go together')
        end
    end  
end

if ~exist('saveGapInfo', 'var')
    saveGapInfo.isGap = 0; 
end
toc(t)


%% save intan info
intaninfo.amplifier_channels = intandata{i}.amplifier_channels;
intaninfo.fileinfo = intandata{i}.fileinfo;
intaninfo.datapts = sum(datapointsrec);%clear datapointsrec;
intaninfo.nAmpChannels = nAmpChannels;
intaninfo.nDigInChannels = nDigInChannels;
intaninfo.nADCChannels = nADCChannels;
save(fullfile(destinationfolder, ['intaninfo', num2str(index(j)), '.mat']), 'intaninfo');
%end moving end to after everything is made (line - 409)
%% load file info to get information about desired start and endtime, etc
%load(fullfile(destinationfolder, ['fileinfo', num2str(index(3)), '.mat']))
starttime = NaN;
endtime = NaN;
raweegsamplerate = params.samprate; 

%% make raweeg data structure shell for saving
%%% defaults
applystarttime = 0; applyendtime = 0; 
raweeg.starttime = 0;
raweeg.endtime = [];

%%% info for structure
samprate =  intaninfo.fileinfo.amplifier_sample_rate;
raweeg.descript = 'raweeg data';
raweeg.index = index;
raweeg.samprate = raweegsamplerate;
raweeg.downsample = samprate/raweegsamplerate;

%%% apply start and end times
if ~isnan(endtime) && isnumeric(endtime)
    endindex = endtime*60*samprate; %endtime is in minutes
    if endindex <= sum(datapointsrec)
        applyendtime = 1;
        raweeg.endtime = endtime;
    else
        warning(['endtime is longer than file length for index: ', ...
            num2str(index), '. not applying endtime'])
    end
end

if starttime>0 && isnumeric(starttime)
    startindex = starttime*samprate; %starttime is in seconds
    applystarttime = 1;
    raweeg.starttime = starttime;
end

%% make raweeg20k or raweeg30k data structure shell
eval([['raweeg' num2str(samprate/1000) 'k'],'.samprate = samprate;'])
eval([['raweeg' num2str(samprate/1000) 'k'],'.downsample = 0;'])
eval([['raweeg' num2str(samprate/1000) 'k'],'.gapInfo = saveGapInfo;'])
eval([['raweeg' num2str(samprate/1000) 'k'],'.descript = ''raweeg data not downsampled. all other props same as raweeg'';'])

%% disp update
disp(['data saving for file ', num2str(index(j)) ])

%% load amplifier data and save
tic
for ch = 1:nAmpChannels
    frewind(ampfid{ch}) %move to the beginning of the file for reading
    chdat = fread(ampfid{ch}, [1, sum(datapointsrec)], 'double');
    %chdat = reshape(chdat, [1, sum(datapointsrec)]); %reshape into column 

    %%% check size is correct
    if ~isequal(size(chdat,2), sum(datapointsrec))
        error('something went wrong with temp file writing or reading')
    end

    %%% check reading is happening properly
    if ch == 1
        if ~isequal(chdat(end-100:end), last100)
            error('something went wrong with temp file reading or writing')
        end
    end

    %%% apply start and end times
    if applystarttime
        chdat = chdat(startindex:end);
    end
    if applyendtime
        chdat = chdat(1:endindex);
    end

    %%% raw eeg 20/30k
    eval([['raweeg' num2str(samprate/1000) 'k'],'.data = chdat;'])

    %%% downsample raweeg
    raweegdownsample = samprate/params.samprate;
    if raweegdownsample > 1 && rem(raweegdownsample,1)==0
        raweeg.data = downsample(chdat, raweegdownsample);
    elseif raweegdownsample == 1
        raweeg.data = chdat;
    elseif raweegdownsample>1 && rem(raweegdownsample,1)~=0 %decimete raweeg to match raweegsamplerate
        raweeg.data = resample(chdat, samprate, raweegsamplerate);
        warning(['resampling raweeg data from ', num2str(samprate), ' to achieve proper sampling rate of ', num2str(raweegsamplerate), '.'])
    else
        error ('raweeg sampling rate is higher than data collected')
    end

    %%% what brain region and save channel does this ind belong to? 
    %currently this only works for 2 probes
    %updated ALP 4/8/21 to work for recording configurations that have
    %uneven # of channels, like a 32 ch and 64 ch probe
    p = 1; saveChI = ch; 
    if ch > 1 && (ch/length(params.savechnum{1})) > 1
        p = 2;
        saveChI = ch - length(params.savechnum{1}); 
    end

    %%% add information about the native port and channel of this data
    raweeg.port = intaninfo.amplifier_channels(ch).port_name;
    raweeg.nativechannel = intaninfo.amplifier_channels(ch).native_channel_name;
    raweeg.channelInd = ch; 
    raweeg.portReg = params.brainReg{p};
    raweeg.savechannel = params.savechnum{p}(saveChI);
    raweeg.gapInfo = saveGapInfo; 


    %%% create the proper directories for the channel
    savechdir = [destinationfolder, params.brainReg{p}, '\', num2str(params.savechnum{p}(saveChI)), '\'];
    if ~isfolder(savechdir)
        mkdir(savechdir);
    end

    %%% save raweeg.mat and raweeg20/30k.mat
    tempraweegfilename = ['raweeg', num2str(index(j))];
    save([savechdir, tempraweegfilename], 'raweeg', '-v7.3')

    tempraweegfilename2 = [['raweeg' num2str(samprate/1000) 'k'], num2str(index(j))]; %filename based on variable samprate
    save([savechdir, tempraweegfilename2], ['raweeg' num2str(samprate/1000) 'k'], '-v7.3')

    %%% close and clear
    fclose(ampfid{ch}); 
    clear chdat
end
toc

%% load other data and save
%helpful stuff from params
ttlsamprate = params.ttlsamprate; 

%digital input data
frewind(digfid) %move to the beginning of the file for reading
digdat = fread(digfid, [nDigInChannels, sum(datapointsrec)], 'double');

%will probably need to add something here to deal with times whne there is
%no adc data
frewind(adcfid)
adcdat = fread(adcfid, [nADCChannels, sum(datapointsrec)], 'double');

%%% make TTL structure shells
for t = 1: length(params.ttlfilename)
    eval([params.ttlfilename{t},'.descript = ''ttl data'';'])
    eval([params.ttlfilename{t},'.index = index;'])
    eval([params.ttlfilename{t},'.samprate = ttlsamprate;']) 
    eval([params.ttlfilename{t},'.downsample = samprate/ttlsamprate;'])
end

%%% apply start and endtimes
if ~isnan(endtime) && isnumeric(endtime) %exclude ? cases
    for t = 1: length(params.ttlfilename)
        eval([params.ttlfilename{t},'.endtime = endtime;']) % is in seconds already
    end

    endindex = endtime*60*samprate; %endtime is in minutes
    if endindex <= size(digdat,2)
        digdat = digdat(:,1:endindex);
        adcdat = adcdat(:,1:endindex);
    else
        warning(['endtime is longer than file length for index: ', num2str(index), '. not applying endtime'])
    end
else
    for t = 1: length(params.ttlfilename)
        eval([params.ttlfilename{t},'.endtime = [];'])
    end
end

%starttime
if starttime>0 && isnumeric(starttime)%exclude ? cases
    for t = 1: length(params.ttlfilename)
        eval([params.ttlfilename{t},'.starttime = starttime;'])
    end

    startindex = starttime*samprate; %starttime is in seconds
    digdat = digdat(:,startindex:end);
    adcdat = adcdat(:,startindex:end);
else
    for t = 1: length(params.ttlfilename)
        eval([params.ttlfilename{t},'.starttime = 0;'])%'{index(1)}{index(2)}{index(3)}.starttime = 0;'
    end
end

%%% downsample and save TTL data
for t = 1:length(params.ttlfilename)
    %downsample ttl
    ttldownsample = samprate/ttlsamprate;
    if ttldownsample > 1 && rem(ttldownsample,1)==0
        eval([params.ttlfilename{t},'.data = downsample(digdat(params.ttlCh(t)+1,:), ttldownsample);'])
    elseif ttldownsample == 1
        eval([params.ttlfilename{t},'.data = digdat(params.ttlCh(t)+1,:);'])
    elseif ttldownsample>1 && rem(ttldownsample,1)~=0
        eval([params.ttlfilename{t},'.data = resample(digdat(params.ttlCh(t)+1,:), samprate, ttlsamprate);'])
        warning(['resampling ttl data from ', num2str(samprate), ' to achieve proper sampling rate of ', num2str(ttlsamprate), '.'])
    else
        error ('ttl sampling rate is higher than data collected')
    end
    eval([params.ttlfilename{t},'.digChannel = params.ttlCh(t);'])

    %save
    tempttlfilename = [params.ttlfilename{t}, num2str(index(j))];
    save([destinationfolder,'/', tempttlfilename], params.ttlfilename{t}, '-v7.3')
end

%% save trigger data 
eval([params.triggerfilename,'.autotrigger = digdat(params.autoTriggerCh+1,:);']) 
eval([params.triggerfilename,'.recordingtrigger = digdat(params.recordingTriggerCh+1,:);']) 
eval([params.triggerfilename,'.samprate = samprate;']) 
eval([params.triggerfilename,'.starttime = raweeg.starttime;']) 
eval([params.triggerfilename,'.endtime =  raweeg.endtime;']) 
save([destinationfolder, '/' strcat(params.triggerfilename, num2str(index(j)))], params.triggerfilename, '-v7.3');


%% load ADC data and save

%%% make ADC structure and save
if isfield(params, 'adc_channels')
    velocity.data1 = adcdat(1,:);
    velocity.data2 = adcdat(2,:);
    velocity.samprate = samprate;
    velocity.info = params.adcInfo;
    velocity.starttime = raweeg.starttime;
    velocity.endtime =  raweeg.endtime;
    save([destinationfolder, '/', strcat(params.velocityfilename, num2str(index(j)))], 'velocity', '-v7.3');
end

% if isfield(params, 'driveLED_channel') %%stim LED driver voltage for Nuri 12/23/21
%     stimLED{index(1)}{index(2)}{index(3)}.data = adcdat(3,:);
%     stimLED{index(1)}{index(2)}{index(3)}.samprate = samprate;
%     stimLED{index(1)}{index(2)}{index(3)}.info = params.ledInfo;
%     stimLED{index(1)}{index(2)}{index(3)}.starttime = raweeg{index(1)}{index(2)}{index(3)}.starttime;
%     stimLED{index(1)}{index(2)}{index(3)}.endtime =  raweeg{index(1)}{index(2)}{index(3)}.endtime;
%     save([destinationfolder, '/', strcat(params.LEDfilename, num2str(index(3)))], 'stimLED', '-v7.3');
% end

fclose(adcfid);
fclose(digfid);
%% delete temp files 
%do I want to do this here? 

for ch = 1:length(cell2mat(params.savechnum))
    chfilename = [tempdir, 'tempdata', num2str(index(j)), '_ch', num2str(ch-1), '.bin'];
    delete(chfilename);
end
clear datapointsrec;
% other data types
delete(digfilename);
delete(adcfilename);

cd(destinationfolder)
rmdir tempfiles
cd('\\ad.gatech.edu\bme\labs\singer\Danielle\code\AnalysisCode\Neuropixels_analyses')
end
end
