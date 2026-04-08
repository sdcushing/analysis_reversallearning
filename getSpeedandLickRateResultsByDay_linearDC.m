function [speedDatAll, lickRateDatAll] = getSpeedandLickRateResultsByDay_linearDC(dirs,uniqSess,params)
%code follows similar structure as getBehaviorROC_JLK.m,
% getFiringRateResultsByDay_linearJLK, and getDecodingResultsByDay_linearJLK

%% Create struct for individual days %%
for id = 1:2%length(params.decoding.decID)%3 not yet operational looks like

    %loop through unique sessions
    for ss = 1:size(uniqSess,1)

        dayDatafname = fullfile([dirs.saveoutputstructs, 'Data\Neural\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' params.decoding.decID{id} 'Data.mat']);

        %loop through environments
        for ee = 1:2%:length(params.environments)
            currEnv = params.environments{ee};

            speedfname = fullfile([dirs.saveoutputstructs, 'Data\Behavior\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' currEnv params.decoding.decID{id} 'Speed.mat']);
            lickRatefname = fullfile([dirs.saveoutputstructs, 'Data\Behavior\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' currEnv params.decoding.decID{id} 'LickRate.mat']);

            if isfile(dayDatafname) %&& ~isfile(speedfname)
                %load this day's data
                data = load(dayDatafname);
                data = data.data;

                %only loop if session exists (included because og sessions
                % have fewer entries than update sessions)
                if ~isempty(data.(currEnv))

                    %initialize structs for this session
                    speedDat = [];
                    lickRateDat = [];

                    %collect speed and lick data
                    if strcmp(params.decoding.decID{id}, 'lap')
                        speedDat = squeeze(data.(currEnv).degBinSpeedSmooth(end,:,:));
                        lickRateDat = squeeze(data.(currEnv).degBinLickRateSmooth(end,:,:));
                    elseif strcmp(params.decoding.decID{id}, 'trial')
                        speedDat = squeeze(data.(currEnv).degBinSpeedSmooth(end,:,:,:));
                        lickRateDat = squeeze(data.(currEnv).degBinLickRateSmooth(end,:,:,:));
                    end

                    %save data for this session
                    dayDataDir = fullfile([dirs.saveoutputstructs, 'Data\Behavior\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2))]);
                    if ~isfolder(dayDataDir)
                        mkdir(dayDataDir)
                    end
                    save(speedfname, 'speedDat', '-v7.3');
                    save(lickRatefname, 'lickRateDat', '-v7.3');

                end%if ~isempty(data.(currEnv))
            end%if isfile(dayDatafname) && ~isfile(lapRateMapfname)

        end%ee

        sprintf('Finished gathering speed and lick rate per %s per day for %s%d_%d', params.decoding.decID{id}, params.iden, uniqSess(ss,1), uniqSess(ss,2))

    end%ss

end%id

%% Create struct for all days %%
%done separately from individual days due to computational restraints
for id = 1:2%:length(params.decoding.decID)

    %initialize output structs
    if strcmp(params.decoding.decID{id}, 'lap')
        %output = sessions x environments x lap x position bin
        speedDatAll = nan(size(uniqSess,1),2,200,180);
        lickRateDatAll = nan(size(uniqSess,1),2,200,180);
    elseif strcmp(params.decoding.decID{id}, 'trial')
        %output = sessions x environments x zone type x trial x position bin
        speedDatAll = nan(size(uniqSess,1),2,3,400,50);
        lickRateDatAll = nan(size(uniqSess,1),2,3,400,50);
    end%if lap/trial

    %loop through unique sessions
    for ss = 1:size(uniqSess,1)

        %loop through environments
        for ee = 1:2%:length(params.environments)
            currEnv = params.environments{ee};

            speedfname = fullfile([dirs.saveoutputstructs, 'Data\Behavior\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' currEnv params.decoding.decID{id} 'Speed.mat']);
            lickRatefname = fullfile([dirs.saveoutputstructs, 'Data\Behavior\dayData\' params.iden num2str(uniqSess(ss,1)) '\' num2str(uniqSess(ss,2)) '\' currEnv params.decoding.decID{id} 'LickRate.mat']);

            if isfile(speedfname)

                %load this day's data
                speedDat = load(speedfname);
                speedDat = speedDat.speedDat;
                lickRateDat = load(lickRatefname);
                lickRateDat = lickRateDat.lickRateDat;

                %collect speed and lick data
                if strcmp(params.decoding.decID{id}, 'lap')
                    %add data to struct for all sessions
                    % output = sessions x environments x lap x position bin
                    speedDatAll(ss,ee,1:size(speedDat,1),1:size(speedDat,2)) = speedDat;
                    lickRateDatAll(ss,ee,1:size(lickRateDat,1),1:size(lickRateDat,2)) = lickRateDat;
                elseif strcmp(params.decoding.decID{id}, 'trial')
                    %add data to struct for all sessions
                    % output = sessions x environments x zone type x trial x position bin
                    speedDatAll(ss,ee,1:size(speedDat,1),1:size(speedDat,2),1:size(speedDat,3)) = speedDat;
                    lickRateDatAll(ss,ee,1:size(lickRateDat,1),1:size(lickRateDat,2),1:size(lickRateDat,3)) = lickRateDat;
                end

            end%if isfile(speedfname)
        end%ee

        sprintf('Finished gathering speed and lick rate per %s all days', params.decoding.decID{id})

    end%ss

    %save data
    speedAllfname = fullfile([dirs.saveoutputstructs, 'Data\Behavior\dayData\' params.decoding.decID{id} 'SpeedAll.mat']);
    lickRateAllfname = fullfile([dirs.saveoutputstructs, 'Data\Behavior\dayData\' params.decoding.decID{id} 'LickRateAll.mat']);
    save(speedAllfname, 'speedDatAll', '-v7.3');
    save(lickRateAllfname, 'lickRateDatAll', '-v7.3');

end%id

end%function
