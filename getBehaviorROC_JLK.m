function [ROC] = getBehaviorROC_JLK(allindex,dirs,uniqSess,params)
%adapted from getNovelBehaviorROC_JLK and getBehaviorDistributionAZvsNRZ_JLK

%% create or load ROC data for each session %%

ROC = [];
sessionInfo = [];

for id = 1:length(params.rocID)
    %loop through unique sessions
    for ss = 1:size(uniqSess,1)

        animal = uniqSess(ss, 1);
        sessDate = uniqSess(ss, 2);
        tempidx = find(ismember(allindex(:,1), animal) & ismember(allindex(:,2), sessDate) & ismember(allindex(:,6), 2));%include this animal, this date, active only
        sessionInfo = allindex(tempidx, :);
        temp{ss,1} = sessionInfo;

        ROCfname = fullfile(dirs.saveoutputstructs, 'Data\Behavior\ROC', [params.iden num2str(animal)], ...
            num2str(sessionInfo(1,2)), [params.rocID{id} '.mat']);

        %get latest file with name that contains the specified string if
        % exist, make speed/lickrate AZ vs NRZ distributiuon file if not
        % adapted from getBehaviorDistributionAZvsNRZ_JLK
        if ~exist(ROCfname) || params.rewrite.ROC

            %%%%% make save path if not already created %%%%%
            if ~isfolder([dirs.saveoutputstructs, 'Data\Behavior\ROC'])
                mkdir([dirs.saveoutputstructs, 'Data\Behavior\ROC'])
            end

            %specify some variables
            if strcmp(params.rocID{id}, 'speed')
                fieldname2select = 'velocCountsSmooth';
            elseif strcmp(params.rocID{id}, 'lickrate') || strcmp(params.rocID{id}, 'deltalickrate')
                fieldname2select = 'lickRateSmooth';
            else
                print('Select from behavior type options: speed, lickrate')
            end
            fNames = {'az','cz', 'nevcz', 'sessionInfo'};

            %initialize data structure with empty fields
            data = [];
            for ee = 1:length(params.environments)
                for fn = 1:length(fNames)
                    data.(params.environments{ee}).(fNames{fn}) = [];
                end
            end

            %loop through files
            files = sessionInfo(:,3);
            for f = 1:length(files)
                currFile = files(f);
                sessionType = sessionInfo(f,4);

                datafname = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData\' params.iden num2str(animal)], ...
                    [num2str(sessDate) '_' num2str(currFile) '_' num2str(sessionType)], 'statsByLap.mat');
                if isfile(datafname)
                    load(datafname)
                end

                if ~isempty(statsByLap)
                    if sessionType == 1 %session in original environment
                        currEnv = 'og';
                    elseif sessionType == 2 %session in update environment
                        currEnv = 'up';
                    elseif sessionType == 3 %session in novel environment
                        currEnv = 'nov';
                    elseif sessionType == 4 %session in novel2 environment
                        currEnv = 'nov2';
                    else
                        disp(['Session type not specified for day ' num2str(sessionInfo(f,2))])
                    end

                    %get trial data to select AZ/CZ from
                    lapData = statsByLap.(fieldname2select);
                    if strcmp(params.rocID{id}, 'deltalickrate')
                        %JLK fixed 4/17/25 to add column of zeros then subtract (previously subtracted then added column of zeros)
                        lapData = [zeros(size(lapData,1),1) lapData];
                        lapData = diff(lapData,1,2);
                    end

                    %get zone info
                    [params.Azones, params.Rzones, params.NRzones, params.NevRzones] = getZoneInfo_linearJLK(statsByLap.fileInfo, [params.iden num2str(animal)]);

                    for lp = 1:size(lapData,1)
                        for zn = 1:length(params.Azones)
                            tmpBins = [];
                            tmpBins = params.Azones(zn):params.binsize_deg:params.Azones(zn)+statsByLap.fileInfo.cueSize-params.binsize_deg;
                            tmpBins = round(tmpBins/params.binsize_deg);%DC added round to handle offset zones with new projector
                            data.(currEnv).az = [data.(currEnv).az; lapData(lp,tmpBins)];
                        end
                        for zn = 1:length(params.NRzones)
                            tmpBins = [];
                            tmpBins = params.NRzones(zn):params.binsize_deg:params.NRzones(zn)+statsByLap.fileInfo.cueSize-params.binsize_deg;
                            tmpBins = round(tmpBins/params.binsize_deg);
                            data.(currEnv).cz = [data.(currEnv).cz; lapData(lp,tmpBins)];
                        end
                        for zn = 1:length(params.NevRzones)
                            if ~isnan(params.NevRzones(zn))
                                tmpBins = [];
                                tmpBins = params.NevRzones(zn):params.binsize_deg:params.NevRzones(zn)+statsByLap.fileInfo.cueSize-params.binsize_deg;
                                tmpBins = round(tmpBins/params.binsize_deg);
                                data.(currEnv).nevcz = [data.(currEnv).nevcz; lapData(lp,tmpBins)];
                            else
                                data.(currEnv).nevcz = [];
                            end
                        end
                    end

                    %save info common to data structure
                    data.(currEnv).azBins_deg = params.Azones;
                    data.(currEnv).czBins_deg = params.NRzones;
                    data.(currEnv).nevczBins_deg = params.NevRzones;
                    data.(currEnv).sessionInfo = [data.(currEnv).sessionInfo; sessionInfo(f,:)];

                end
            end

            %save data
            dir2save = fullfile(dirs.saveoutputstructs, 'Data\Behavior\ROC', [params.iden num2str(animal)], ...
                num2str(sessionInfo(1,2)));
            if ~isfolder(dir2save); mkdir(dir2save); end
            save([dir2save, '\', params.rocID{id}, '.mat'], 'data', 'data');


        else
            load(ROCfname);
        end


        %% ROC type 1: use all trials and anticipatory zones vs primary control zones %%
        for ee = 1:length(params.environments)
            currEnv = params.environments{ee};
            if ~isempty(data.(currEnv).az) %if session exists
                %reshape to Nx1 structure and combine AZ and CZ data into one
                data_az = nanmean(data.(currEnv).az, 2);
                data_cz = nanmean(data.(currEnv).cz, 2);
                rocData = calcBehaviorROC(data_az*params.rocMultiplier(id), data_cz*params.rocMultiplier(id));

                rocOut.(currEnv).mdl{ss} = rocData.mdl;
                rocOut.(currEnv).scores{ss} = rocData.scores;
                rocOut.(currEnv).X{ss} = rocData.X;
                rocOut.(currEnv).Y{ss} = rocData.Y;
                rocOut.(currEnv).T{ss} = rocData.T;
                rocOut.(currEnv).AUC(ss) = rocData.AUC;
            else
                rocOut.(currEnv).AUC(ss) = nan;
            end
        end

        %% ROC type 2: separate first-half and second-half of trials with anticipatory zones vs primary control zones
        for ee = 1:length(params.environments)
            currEnv = params.environments{ee};
            if ~isempty(data.(currEnv).az) %if session exists
                nTrialsAZ = size(data.(currEnv).az,1); half_az = ceil(nTrialsAZ/2);
                nTrialsCZ = size(data.(currEnv).cz,1); half_cz = ceil(nTrialsCZ/2);
                nTrialsnevCZ = size(data.(currEnv).nevcz,1); half_nevcz = ceil(nTrialsnevCZ/2);

                for iH = 1:2 %first or second half
                    if iH == 1
                        data_az = nanmean(data.(currEnv).az(1:half_az,:), 2);
                        data_cz = nanmean(data.(currEnv).cz(1:half_cz,:), 2);
                        data_nevcz = nanmean(data.(currEnv).nevcz(1:half_nevcz,:), 2);
                    else
                        data_az = nanmean(data.(currEnv).az(half_az+1:end,:), 2);
                        data_cz = nanmean(data.(currEnv).cz(half_cz+1:end,:), 2);
                        data_nevcz = nanmean(data.(currEnv).nevcz(half_cz+1:end,:), 2);
                    end
                    %anticipatory zones vs primary control zones
                    rocData = calcBehaviorROC(data_az*params.rocMultiplier(id), data_cz*params.rocMultiplier(id));
                    rocOut.([currEnv, '_', num2str(iH)]).mdl{ss} = rocData.mdl;
                    rocOut.([currEnv, '_', num2str(iH)]).scores{ss} = rocData.scores;
                    rocOut.([currEnv, '_', num2str(iH)]).X{ss} = rocData.X;
                    rocOut.([currEnv, '_', num2str(iH)]).Y{ss} = rocData.Y;
                    rocOut.([currEnv, '_', num2str(iH)]).T{ss} = rocData.T;
                    rocOut.([currEnv, '_', num2str(iH)]).AUC(ss) = rocData.AUC;
                    if strcmp(params.environments{ee},'up')
                        %update anticipatory zones vs nevernative control zones
                        rocData = calcBehaviorROC(data_az*params.rocMultiplier(id), data_nevcz*params.rocMultiplier(id));
                        rocOut.([currEnv, 'nev_', num2str(iH)]).mdl{ss} = rocData.mdl;
                        rocOut.([currEnv, 'nev_', num2str(iH)]).scores{ss} = rocData.scores;
                        rocOut.([currEnv, 'nev_', num2str(iH)]).X{ss} = rocData.X;
                        rocOut.([currEnv, 'nev_', num2str(iH)]).Y{ss} = rocData.Y;
                        rocOut.([currEnv, 'nev_', num2str(iH)]).T{ss} = rocData.T;
                        rocOut.([currEnv, 'nev_', num2str(iH)]).AUC(ss) = rocData.AUC;
                        %original anticipatory zones vs nevernative control zones
                        rocData = calcBehaviorROC(data_cz*params.rocMultiplier(id), data_nevcz*params.rocMultiplier(id));
                        rocOut.([currEnv, 'nev2_', num2str(iH)]).mdl{ss} = rocData.mdl;
                        rocOut.([currEnv, 'nev2_', num2str(iH)]).scores{ss} = rocData.scores;
                        rocOut.([currEnv, 'nev2_', num2str(iH)]).X{ss} = rocData.X;
                        rocOut.([currEnv, 'nev2_', num2str(iH)]).Y{ss} = rocData.Y;
                        rocOut.([currEnv, 'nev2_', num2str(iH)]).T{ss} = rocData.T;
                        rocOut.([currEnv, 'nev2_', num2str(iH)]).AUC(ss) = rocData.AUC;
                    end%if strcmp(params.environments{ee},'up')

                end%iH

            else
                for iH = 1:2
                    rocOut.([currEnv, '_', num2str(iH)]).AUC(ss) = nan;
                    rocOut.([currEnv, 'nev_', num2str(iH)]).AUC(ss) = nan;
                    rocOut.([currEnv, 'nev2_', num2str(iH)]).AUC(ss) = nan;
                end
            end
        end

        %% ROC type 3: separate 1st + 2nd block of trials with anticipatory zones vs atlernative control zones and primary control zones vs. nevernative control zones
        for ee = 2 %only update sessions
            currEnv = params.environments{ee};
            if ~isempty(data.(currEnv).az) && size(data.(currEnv).az,1) %if session exists and has at least 10 trials
                %get trials to use
                nTrialsAZ = size(data.(currEnv).az,1);
                nTrialsCZ = size(data.(currEnv).cz,1);
                nTrialsnevCZ = size(data.(currEnv).nevcz,1);
                for bl = 1:2
                    %initialize
                    useAzTrials = nan(1,params.numTrPerBlock);
                    useCzTrials = nan(1,params.numTrPerBlock);
                    usenevCzTrials = nan(1,params.numTrPerBlock);
                    if bl == 1
                        useAzTrials = [1:params.numTrPerBlock];
                        useCzTrials = [1:params.numTrPerBlock];
                        usenevCzTrials = [1:params.numTrPerBlock];
                    elseif bl == 2
                        useAzTrials = [nTrialsAZ-params.numTrPerBlock+1:nTrialsAZ];
                        useCzTrials = [nTrialsCZ-params.numTrPerBlock+1:nTrialsCZ];
                        usenevCzTrials = [nTrialsnevCZ-params.numTrPerBlock+1:nTrialsnevCZ];
                    end%if bl == 1

                    %get data
                    data_az = nanmean(data.(currEnv).az(useAzTrials,:), 2);
                    data_cz = nanmean(data.(currEnv).cz(useCzTrials,:), 2);
                    data_nevcz = nanmean(data.(currEnv).nevcz(usenevCzTrials,:), 2);
                    %anticipatory zones vs original anticipatory zones
                    rocData = calcBehaviorROC(data_az*params.rocMultiplier(id), data_cz*params.rocMultiplier(id));
                    rocOut.(sprintf('%s_block%d', currEnv,bl)).mdl{ss} = rocData.mdl;
                    rocOut.(sprintf('%s_block%d', currEnv,bl)).scores{ss} = rocData.scores;
                    rocOut.(sprintf('%s_block%d', currEnv,bl)).X{ss} = rocData.X;
                    rocOut.(sprintf('%s_block%d', currEnv,bl)).Y{ss} = rocData.Y;
                    rocOut.(sprintf('%s_block%d', currEnv,bl)).T{ss} = rocData.T;
                    rocOut.(sprintf('%s_block%d', currEnv,bl)).AUC(ss) = rocData.AUC;
                    %anticipatory zones vs nevernative (never rewarded) control zones
                    rocData = calcBehaviorROC(data_az*params.rocMultiplier(id), data_nevcz*params.rocMultiplier(id));
                    rocOut.(sprintf('%snev_block%d', currEnv,bl)).mdl{ss} = rocData.mdl;
                    rocOut.(sprintf('%snev_block%d', currEnv,bl)).scores{ss} = rocData.scores;
                    rocOut.(sprintf('%snev_block%d', currEnv,bl)).X{ss} = rocData.X;
                    rocOut.(sprintf('%snev_block%d', currEnv,bl)).Y{ss} = rocData.Y;
                    rocOut.(sprintf('%snev_block%d', currEnv,bl)).T{ss} = rocData.T;
                    rocOut.(sprintf('%snev_block%d', currEnv,bl)).AUC(ss) = rocData.AUC;
                    %original reward zones vs nevernative (never rewarded) control zones
                    rocData = calcBehaviorROC(data_cz*params.rocMultiplier(id), data_nevcz*params.rocMultiplier(id));
                    rocOut.(sprintf('%snev2_block%d', currEnv,bl)).mdl{ss} = rocData.mdl;
                    rocOut.(sprintf('%snev2_block%d', currEnv,bl)).scores{ss} = rocData.scores;
                    rocOut.(sprintf('%snev2_block%d', currEnv,bl)).X{ss} = rocData.X;
                    rocOut.(sprintf('%snev2_block%d', currEnv,bl)).Y{ss} = rocData.Y;
                    rocOut.(sprintf('%snev2_block%d', currEnv,bl)).T{ss} = rocData.T;
                    rocOut.(sprintf('%snev2_block%d', currEnv,bl)).AUC(ss) = rocData.AUC;
                end

            else
                for bl = 1:2
                    rocOut.(sprintf('%s_block%d', currEnv,bl)).AUC(ss) = nan;
                    rocOut.(sprintf('%snev_block%d', currEnv,bl)).AUC(ss) = nan;
                    rocOut.(sprintf('%snev2_block%d', currEnv,bl)).AUC(ss) = nan;
                end
            end
        end


        %% ROC type 4: use all trials and anticipatory zones vs nevernative control zones
        for ee = 1:length(params.environments)
            currEnv = params.environments{ee};
            if isfield(data.(currEnv), 'nevcz') && ~isempty(data.(currEnv).nevcz)
                %reshape to Nx1 structure and combine AZ and CZ data into one
                data_az = nanmean(data.(currEnv).az, 2);
                data_nevcz = nanmean(data.(currEnv).nevcz, 2);
                rocData = calcBehaviorROC(data_az*params.rocMultiplier(id), data_nevcz*params.rocMultiplier(id));

                rocOut.(sprintf('%snev', currEnv)).mdl{ss} = rocData.mdl;
                rocOut.(sprintf('%snev', currEnv)).scores{ss} = rocData.scores;
                rocOut.(sprintf('%snev', currEnv)).X{ss} = rocData.X;
                rocOut.(sprintf('%snev', currEnv)).Y{ss} = rocData.Y;
                rocOut.(sprintf('%snev', currEnv)).T{ss} = rocData.T;
                rocOut.(sprintf('%snev', currEnv)).AUC(ss) = rocData.AUC;
            else
                rocOut.(sprintf('%snev', currEnv)).AUC(ss) = nan;
            end
        end


        %% ROC type 5: use all trials and control zones vs nevernative control zones
        for ee = 1:length(params.environments)
            currEnv = params.environments{ee};
            if  isfield(data.(currEnv), 'nevcz') && ~isempty(data.(currEnv).nevcz)
                %reshape to Nx1 structure and combine AZ and CZ data into one
                data_cz = nanmean(data.(currEnv).cz, 2);
                data_nevcz = nanmean(data.(currEnv).nevcz, 2);
                rocData = calcBehaviorROC(data_cz*params.rocMultiplier(id), data_nevcz*params.rocMultiplier(id));

                rocOut.(sprintf('%snev2', currEnv)).mdl{ss} = rocData.mdl;
                rocOut.(sprintf('%snev2', currEnv)).scores{ss} = rocData.scores;
                rocOut.(sprintf('%snev2', currEnv)).X{ss} = rocData.X;
                rocOut.(sprintf('%snev2', currEnv)).Y{ss} = rocData.Y;
                rocOut.(sprintf('%snev2', currEnv)).T{ss} = rocData.T;
                rocOut.(sprintf('%snev2', currEnv)).AUC(ss) = rocData.AUC;
            else
                rocOut.(sprintf('%snev2', currEnv)).AUC(ss) = nan;
            end
        end

        %% Save all data into a giant structure
        ROC.sessInfo = temp;
        ROC.og_all = rocOut.og;
        ROC.og_1stHalf = rocOut.og_1;
        ROC.og_2ndHalf = rocOut.og_2;
        ROC.up_all = rocOut.up;
        ROC.up_1stHalf = rocOut.up_1;
        ROC.up_2ndHalf = rocOut.up_2;
        ROC.up_block1 = rocOut.up_block1;
        ROC.up_block2 = rocOut.up_block2;
        ROC.upnev_all = rocOut.upnev;
        ROC.upnev_1stHalf = rocOut.upnev_1;
        ROC.upnev_2ndHalf = rocOut.upnev_2;
        ROC.upnev_block1 = rocOut.upnev_block1;
        ROC.upnev_block2 = rocOut.upnev_block2;
        ROC.upnev2_all = rocOut.upnev2;
        ROC.upnev2_1stHalf = rocOut.upnev2_1;
        ROC.upnev2_2ndHalf = rocOut.upnev2_2;
        ROC.upnev2_block1 = rocOut.upnev2_block1;
        ROC.upnev2_block2 = rocOut.upnev2_block2;
        ROC.nov_all = rocOut.nov;
        ROC.nov_1stHalf = rocOut.nov_1;
        ROC.nov_2ndHalf = rocOut.nov_2;
        ROC.nov2_all = rocOut.nov2;
        ROC.nov2_1stHalf = rocOut.nov2_1;
        ROC.nov2_2ndHalf = rocOut.nov2_2;

        filedir = fullfile(dirs.saveoutputstructs, 'Data\Behavior\ROC');

        %save output structure
        savefname = fullfile(filedir, ['ROC' '_' params.rocID{id} '_' datestr(now,'yymmdd')]);
        save(savefname, 'ROC', 'uniqSess');


    end%ss
end%id

end%function
