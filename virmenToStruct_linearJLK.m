function [rawDataBySession, virmen_fileInfo] = virmenToStruct_linearJLK(virmenDataPath, saveBehaviorPath)

%load virmen data
virmen_rawData = load(virmenDataPath);
virmen_rawData = virmen_rawData.saveDatCopy;

%find columns to create data structure
timeColumn = find(strcmp(virmen_rawData.dataHeaders,'time'));
thetaColumn = find(strcmp(virmen_rawData.dataHeaders,'theta'));
rewardsColumn = find(strcmp(virmen_rawData.dataHeaders,'numRewards'));
licksColumn = find(strcmp(virmen_rawData.dataHeaders,'numLicks'));
zoneColumn = find(strcmp(virmen_rawData.dataHeaders,'currentZone'));
transVeloColumn = find(strcmp(virmen_rawData.dataHeaders,'transVeloc'));
rotVeloColumn = find(strcmp(virmen_rawData.dataHeaders,'rotVeloc'));
ttlColumn = find(strcmp(virmen_rawData.dataHeaders,'AutoTrigger'));

%create data structure
rawDataBySession = struct('vrTime',(virmen_rawData.data(:,timeColumn) - virmen_rawData.data(1,timeColumn))*24*3600,...
    'currentDeg',virmen_rawData.data(:,thetaColumn),'rewards',virmen_rawData.data(:,rewardsColumn),...
    'licks',virmen_rawData.data(:,licksColumn),'currentZone',virmen_rawData.data(:,zoneColumn), ...
    'transVelo', virmen_rawData.data(:,transVeloColumn), 'rotVelo', virmen_rawData.data(:,rotVeloColumn), ...
    'ttl', virmen_rawData.data(:,ttlColumn));

%save basic file info
virmen_fileInfo = rmfield(virmen_rawData,{'data'});
if isfield(virmen_fileInfo, 'datalicks')
    virmen_fileInfo = rmfield(virmen_fileInfo, 'datalicks');
end
if isfield(virmen_fileInfo, 'dataHeaders')
    virmen_fileInfo = rmfield(virmen_fileInfo, 'dataHeaders');
end
if ~isfield(virmen_fileInfo, 'cueSize')
    virmen_fileInfo.cueSize = 10; %length of a zone in degrees
end
if isfield(virmen_fileInfo, 'thetaReward')
    virmen_fileInfo.degReward = virmen_fileInfo.thetaReward;
    virmen_fileInfo = rmfield(virmen_fileInfo, 'thetaReward');
end

%% Save data %%
filename = [saveBehaviorPath '\' 'rawDataBySession.mat'];
save(filename, 'rawDataBySession');

end%function