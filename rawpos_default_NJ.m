function [] = rawpos_default_NJ(index, vrDataFolder, anprocesseddatadir, params)
%% Raw Position Data Structure Construction
%Default rawpos extraction function used by NJ. Integrated for both Intan
%and SpikeGadgets acquisition systems
%NJ 10/12/2021
%NJ edit 3/4/22 to include a new field "stimVoltage" indicating analog
%voltage signal to drive the LEDs and exclude the redundant field
%"rewardPump"

%% ---------- Import VR Data and Array of Synced Triggers ------------
disp(['extracting raw position index: ', num2str(index)])
%import SaveDatCopy from VR file
if isfile(fullfile(vrDataFolder, 'dataWithLickometer.mat'))
    load(fullfile(vrDataFolder, 'dataWithLickometer.mat'));
else
    return
end

%gets data matrix from struct, if saveDatCopy was saved as struct
if isa(saveDatCopy,'struct')
    virmenStruct = saveDatCopy.data;
end

%use data headers to find virmen column info
if isfield(saveDatCopy, 'dataHeaders')
    dataHeaders = saveDatCopy.dataHeaders;
    vrTimes(:,1) = virmenStruct(:, find(strcmp(dataHeaders, 'time')) );
else
    vrTimes(:,1) = virmenStruct(:,1);
end

%converts VR computer times to time elapsed in second since VR start
vrTimes(:,2) = (vrTimes(:,1) - vrTimes(1,1)) * 24 * 60 * 60;


%returns array of synced triggers regardless of acquisition system
syncedArray = triggerNonperiodicSync(index, vrDataFolder, anprocesseddatadir);


%% ------------ Create & Save Raw Position Data Structure -------------
syncStart = find(~isnan(syncedArray(:,2)));
iStart = syncStart(1); %virmen index of the first sync'ed time
iEnd = syncStart(end);

%rawpos structure has all relevant virmen data from the first trig pulse to
%the end - NJ 11.04.19
if isfield(saveDatCopy, 'dataHeaders')
    
    rawpos{index(1)}{index(2)}{index(3)} = struct(...
        'vrCompTimes',vrTimes(iStart:iEnd,1),...
        'vrInd', syncStart,...
        'vrTime',vrTimes(iStart:iEnd,2)-vrTimes(iStart,2),...
        'ephysInd', syncedArray(iStart:iEnd,2),...
        'ephysTimes',syncedArray(iStart:iEnd,3),...
        'currentTheta',virmenStruct(iStart:iEnd,find(strcmp(dataHeaders, 'theta'))),...   
        'rotVel',virmenStruct(iStart:iEnd, find(strcmp(dataHeaders, 'rotVeloc'))),...
        'transVel',virmenStruct(iStart:iEnd, find(strcmp(dataHeaders, 'transVeloc'))),...
        'licks',virmenStruct(iStart:iEnd,find(strcmp(dataHeaders, 'numLicks'))),...
        'rewards',virmenStruct(iStart:iEnd,find(strcmp(dataHeaders, 'numRewards'))),...
        'currentZone',virmenStruct(iStart:iEnd,find(strcmp(dataHeaders, 'currentZone'))),...
        'stimVoltage',virmenStruct(iStart:iEnd,find(strcmp(dataHeaders, 'stimVoltage'))));
else %use old NJ column structure prior to 10/2021
    rawpos{index(1)}{index(2)}{index(3)} = struct(...
        'vrCompTimes',vrTimes(iStart:iEnd,1),...
        'vrInd', syncStart,...
        'vrTime',vrTimes(iStart:iEnd,2)-vrTimes(iStart,2),...
        'ephysInd', syncedArray(iStart:iEnd,2),...
        'ephysTimes',syncedArray(iStart:iEnd,3),...
        'currentTheta',virmenStruct(iStart:iEnd,2),...
        'rotVel',virmenStruct(iStart:iEnd,3),...
        'vel',virmenStruct(iStart:iEnd,4),...
        'licks',virmenStruct(iStart:iEnd,6),...
        'rewards',virmenStruct(iStart:iEnd,5),...
        'currentZone',virmenStruct(iStart:iEnd,7));
end

%this is a sanity check for anyone who ever needs one SP 11.13.18
%figure; plot(trigList,syncedArray(:,2),'o'); hold on;
%plot((rawpos{index(1)}{index(2)}{index(3)}.vrtimes-rawpos{index(1)}{index(2)}{index(3)}.vrtimes(1))*60*60*24,rawpos{index(1)}{index(2)}{index(3)}.indices)

save([anprocesseddatadir,'rawpos',num2str(index(3))],'rawpos');
end
