function syncedArray = triggerNonperiodicSync(index, vrFolder, datadir)
%% triggerSync: syncs Intan trigger times to VR trigger times - MKA
%   MUST RUN FUNCTION AS ADMINISTRATOR in MATLAB 2015
%   Outputs an array listing the synced time/samples that each file records
%   the triggers turning OFF or ON
%   Both intanFile and saveDatCopy should be either loaded in the workspace or
%   in path.
%   intanFile may need to be extracted using read_Intan_RHD2000_file.m
%   The name of the intan file may look like "board_dig_in_data" or tempdata.
%   The name of the matlab file may look like "saveDatCopy".

% NJ default function for synchrnonizing virmen to ephys regardless of
% acqusition system used. Nonperiodic sync pulses were used to match VR
% times and ephys times. 10/12/2021

%% Import trigger data from SpikeGadgets/Intan
trigFile = fullfile(datadir, ['sync' num2str(index(3)) '.mat']);
load(trigFile)

trigger = sync{index(1)}{index(2)}{index(3)};
if isfield(trigger, 'samprate')
    samprate = trigger.samprate;
else
    samprate = 30000; %default is 30000 Hz for Nuri's data
end

if isfield(trigger, 'autotrigger') %Intan system uses digital sync pulses
    intanAutoTrigs = trigger.autotrigger;
    pulseOnIdx = find(intanAutoTrigs)';
elseif isfield(trigger, 'data') %SpikeGadgets system uses analog sync pulses
    sgTime = double(trigger.time); %SpikeGadgets system time for each data point
    %multiply by voltage scaling factor to get uV, then divide by 1000 to get V
    SGautoTrigs = double(trigger.data) .* (0.1950 / 1000);
    
    %if input voltage is greater than 3 V, then print it as 1, and 0 otherwise
    SGautoTrigs(SGautoTrigs>3) = 1;
    SGautoTrigs(SGautoTrigs~=1) = 0;
    pulseOnIdx = find(SGautoTrigs);
end

%% Import behavioral data structure from Virmen
load(fullfile(vrFolder, 'dataWithLickometer'));
vrData = saveDatCopy.data;
%import nonperiodic trigger data and time based on column header if exists
if isfield(saveDatCopy, 'dataHeaders')
    dataHeaders = saveDatCopy.dataHeaders;
    vrTimes = vrData(:, find(strcmp(dataHeaders, 'time')) ) *24*60*60;
    vrAutotrigs = vrData(:, find(strcmp(dataHeaders, 'AutoTrigger')) );
    vrTrigOnOff = vrData(:, find(strcmp(dataHeaders, 'StartEndTrigger')) );
    vrTrigOnOff = find(vrTrigOnOff);
else
    vrTimes = (vrData(:,1)-vrData(1,1))*24*60*60;
    vrAutotrigs = vrData(:,10);
    vrTrigOnOff = find(vrData(:,9));
end
vrTimes = vrTimes - vrTimes(1);
if sum(diff(vrTrigOnOff)>2) > 0
    disp('This recording has more than one trigger on/off period');
end

%% find indices of pulse on times for both virmen and ephys
%indices for autotrigger ON's stored in VR to match
vr = find(vrAutotrigs);
vr2 = vr((diff(vr)>1)); %this step is prob unnecessary since I make the two stim pulses to be at least 3 virmen iterations apart
vr2 = [vr(1); vr2; vr(end)];
vr2 = unique(vr2);

%bc ephys computer marks a lot more 1's to represent a single pulse, only
%find distinct pulses separated by enough samples
ephys2 = pulseOnIdx(diff(pulseOnIdx) > 300);
ephys2 = [pulseOnIdx(1); ephys2; pulseOnIdx(end)];
ephys2 = unique(ephys2);

%align virmen and ephys time gap signals to find the best match
timegap_vr = diff(vrTimes(vr2));
timegap_ephys = diff(ephys2 ./ 30000);
nDelay = finddelay(timegap_vr, timegap_ephys, 5); %determine the number of samples one signal precedes/follows
if nDelay > 0 %first signal leads
    ephys2 = ephys2(nDelay+1:end);
    timegap_ephys = timegap_ephys(nDelay+1:end);
elseif nDelay < 0 %second signal leads
    vr2 = vr2(-nDelay+1:end);
    timegap_vr = timegap_vr(-nDelay+1:end);
end

%in rare occasion ephys pulse length is on for much longer than virmen
%detected - in this case, use whatever is shorter and cut off the overhang
if length(ephys2) ~= length(vr2)
    disp('Lengths do not match. Applying minimum length to virmen and ephys.')
    minLength = min([length(ephys2), length(vr2)]);
    ephys2 = ephys2(1:minLength);
    vr2 = vr2(1:minLength);
    
    minLength = min([length(timegap_ephys), length(timegap_vr)]);
    timegap_ephys = timegap_ephys(1:minLength);
    timegap_vr = timegap_vr(1:minLength);
end

%ensure pulse time gap difference between virmen and ephys is reasonable,
%ideally, time difference should be a very small number
pulse_incl = abs(timegap_vr - timegap_ephys) < 0.05;
ephys2 = ephys2([true; pulse_incl]);
vr2 = vr2([true; pulse_incl]);

%report time difference in seconds between the start and end of autotrigger
%detected by VR and Ephys computers
start2end_VR = vrTimes(vr2(end)) - vrTimes(vr2(1));
start2end_Ephys = (ephys2(end) - ephys2(1)) ./ samprate;
disp(['Difference in autotrigger duration between Virmen and Ephys: ' ...
    num2str(abs( start2end_Ephys - start2end_VR )) ' sec.'])


%% Synchronize each Virmen sample with SpikeGadgets data point

%initialize synchronized data structure - synchronization starts only after
%the very first trigger pulse initated by the key press (F3)
syncedArray = nan(size(vrData,1), 3);

%first column is virmen computer time to sync
syncedArray(:,1) = saveDatCopy.data(:,1);

%second column is for indices of ephys data
for i = 1:length(vr2)
    syncedArray(vr2(i),2) = ephys2(i);
end

%test that numbers go up as expecte and report if not
temp = syncedArray(:,2);
temp = temp(~isnan(temp)); %remove NaN's
if sum(temp ~= sort(temp,'ascend')) > 0 
    disp('Ephys-Virman matched values are not in the increasing order.')
end

%time difference from one VR iteration to next; represented in the number of samples by ephys computer
diffVR = [0; diff(syncedArray(:,1)) .* 24 * 60 * 60];
diffVR = round(diffVR .* samprate);

%fill in the ephys index gaps with VR intervals between iterations
for i = 1:length(vr2)-1
    gapNum = minus(vr2(i+1),vr2(i))-1;
    for p = 0:gapNum-1
        prev = syncedArray(vr2(i) + p , 2);
        stepsize = diffVR(vr2(i) + p);
        syncedArray(vr2(i) + p + 1, 2) = prev + stepsize;
    end
end

%third column is SpikeGadgets time if SG system
if exist('sgTime','var')
    syncedArray(vr2(1):vr2(end),3) = sgTime(syncedArray(vr2(1):vr2(end),2));
end


%% plotting 
figure; ax = axes('NextPlot','add','Box','off');
plot(ax, timegap_vr(pulse_incl) );
plot(ax, timegap_ephys(pulse_incl) );
title(num2str(index)); legend('Virmen','Ephys'); ylabel('Time gap between pulses (s)'); xlabel('Pulse index')
