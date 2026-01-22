%created by JLK on 6/10/25
%last check JLK on 6/30/25

%% Set data directories, params, and sessions to use %%
%Note: selectindextable_JLK filter options: animal, id, datesincluded, datesexcluded, SessionNum,
% SessionType (1 = passive, 2 = active, 3 = rest), Track (1 = "TrackA"
% [original], 2 = "TrackA' [update], 3 = "TrackB" [novel], 4 = "TrackC" [novel2]), Recording
%MAJOR NOTE: If using spreadsheet with "VR" column for track info and
% "Recording" column for session number, do the following: 1) Change the
% "VR" column header to "SessionNum" and "VR to "Track"; and 2) copy and paste 
% =IF(G3:G96=1,"TrackA", IF(G3:G96=2,"TrackB", IF(G3:G96=3,"TrackC"))) in a
% cell below the final row then copy the values into the "Track" column to
% replace track numbers with track names (YOU WILL NEED TO CHANGE THE
% G3:G96 VALUES TO ALIGN WITH YOUR ROWS/COLUMN). 
addpath('\\ad.gatech.edu\bme\labs\singer\Danielle\code\vr_novelty_behavior-xiao\commonfunc')

%for Josh's virmen output spreadsheet
[dirs, params] = getDirectoriesAndParams_JLK_DC();%DC update this to call in commonfunc one with User so I only need 1
allindexT = selectindextable_JLK(dirs.spreadsheet, 'animal', params.animals, 'Track', [1 2 3 4]);%'SessionType', [2 3],

% % for Danielle/Xiao's virmen output spreadsheet
% [dirs, params] = getDirectoriesAndParams_JLK_DC();
% allindexT = selectindextable_JLK(dirs.spreadsheet, 'animal', params.animals, 'Track', [1 2 3 4]);

%create allindex and uniqSess variables
allindex = allindexT{:,{'Animal', 'Date', 'SessionNum', 'Track'}};
%add recording info if in spreadsheet, otherwise assume not recording
if sum(strcmp(allindexT.Properties.VariableNames, 'Recording')) ~= 0
    allindex(:,5) = allindexT{:,{'Recording'}};
else
    allindex(:,5) = 0;
end
%add session type info if in spreadsheet, otherwise assume active session
if ~isempty(allindexT{1,'SessionType'})
    allindex(:,6) = allindexT{:,{'SessionType'}};
else
    allindex(:,6) = 2;
end
[uniqSess, ind] = unique(allindex(:,1:2), 'rows'); %define session as one date

%% Specify what you want to analyze here %%
createBehaviorStructs = 0;
plotBehavior = 0;
gatherNeuralData = 1;
doDecoding = 1;

%% Create behavior data structs %%
if createBehaviorStructs

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Lap and Trial Strcuts %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:size(allindex,1) %loop through every session
        addpath('\\ad.gatech.edu\bme\labs\singer\Danielle\code\vr_novelty_behavior\commonfunc')
        addpath('\\ad.gatech.edu\bme\labs\singer\Danielle\code\vr_novelty_behavior\functions')
        %%%%% session info %%%%%
        subj = [params.iden num2str(allindex(i,1))];
        sessDate = num2str(allindex(i,2));
        sessNum = num2str(allindex(i,3));
        processedDataPath = fullfile(dirs.processeddata,[subj '_' sessDate '\']);
        neuralRawDataPath = fullfile(dirs.rawdata, [subj '_' sessDate]);
        trackInfo = num2str(allindex(i,4));
        virmenDataPath = fullfile(dirs.virmenrawdata, ['\' subj '_' sessDate '_' sessNum '\' 'dataWithLickometer.mat']);
        saveBehaviorPath = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData' '\' subj '\' sessDate '_' sessNum '_' trackInfo]);

        if isfile(virmenDataPath)

            %%%%% make save path if not already created %%%%%
            if ~isfolder(saveBehaviorPath)
                mkdir(saveBehaviorPath)
            end

            %%%%% create virmen data structure, file information, and define zones %%%%%%
            %adapted from virmenToStruct_linearNJ.m
            %Note: If using update track, user may need to change the reward and non-reward locations.
            [~, virmen_fileInfo] = virmenToStruct_linearJLK(virmenDataPath, saveBehaviorPath);
            if virmen_fileInfo.sessioninfo == 'DC21_251118_5'
                virmen_fileInfo.trackname = 'TrackA''';
            end
            %rhd2mat_tempbin_DC(neuralRawDataPath, processedDataPath, sessNum, params);
            disp(['Extracting Virmen Data: ', subj, ' ', sessDate, ' ', sessNum])
            anvrdatafolder = fullfile(dirs.virmenrawdata, [subj '_', sessDate, '_',  sessNum]);
            Args = {sessNum, anvrdatafolder, processedDataPath, params};
            %feval(params.exportbehaviorfunc, Args{:});
            rawposfile = fullfile(processedDataPath, sprintf('rawpos%s.mat', sessNum));
            rawpos = load(rawposfile);
            rawDataBySession = rawpos.rawpos;
            save(fullfile(saveBehaviorPath, 'rawDataBySession.mat'), 'rawDataBySession');
            [params.Azones, params.Rzones, params.NRzones, params.NevRzones] = getZoneInfo_linearJLK(virmen_fileInfo, subj);

            %%%%% create behavior structs %%%%%
            if allindex(i,6) == 2%Active VR session

                %LAPS
                %adapted from getTrialByTrialStats_linearJLK and getSessionStats_linearJLK
                %Note: only completed laps
                if ~isfile([saveBehaviorPath '\' 'statsByLap.mat']) || params.rewrite.behavior
                    getLapBehaviorStats_linearDC(rawDataBySession, virmen_fileInfo, params, saveBehaviorPath);
                end

                %TRIALS
                %adapted from getTrialByTrialStats_linearJLK and getSessionStats_linearJLK
                %Note: uses all trials
                if ~isfile([saveBehaviorPath '\' 'statsByRewardTrial.mat']) || params.rewrite.behavior
                    getTrialBehaviorStats_linearDC(rawDataBySession, virmen_fileInfo, params, saveBehaviorPath);
                end

            elseif allindex(i,6) == 3%Rest session
                
                %REST
                %adapted from getRestSessionStats_linearJLK
                if ~isfile([saveBehaviorPath '\' 'statsByRestSession.mat']) || params.rewrite.behavior
                    getRestBehaviorStats_linearJLK(rawDataBySession, virmen_fileInfo, saveBehaviorPath);
                end

            end%if allindex(i,6)

        else
            sprintf('File not found - Skipped analyzing session data for %s_%s_%s', subj, sessDate, sessNum)
        end%isfile(virmenDataPath)

        fprintf(['Finished analyzing trial data for session ' num2str(i) ' out of ' num2str(size(allindex,1)) ': ' params.iden num2str(allindex(i,1)) '_' num2str(allindex(i,2)) '_' num2str(allindex(i,3)) '... \n'])

    end%i

    %%%%%%%%%%%%%%%%
    %%%%% ROCs %%%%%
    %%%%%%%%%%%%%%%%
    addpath('\\ad.gatech.edu\bme\labs\singer\Danielle\code\vr_novelty_behavior\functions')
    
    getBehaviorROC_DC(allindex,dirs,uniqSess,params)

end%if createBehaviorStructs

%% Plot behavior %%
if plotBehavior

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Lap and Trial Strcuts %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %CHOICES:
    % Either 1) plot by date by setting plotByDateOrType equal to 1 and
    %  setting plotDate equal to the date you want to plot (indexed by mouse);
    %  or 2) plot by original, update, novel, or novel2 day by setting plotByDateOrType
    %  equal to 1, 2, 3, or 4, respectively, and setting whichDay equal to the
    %  day number you want to plot (e.g. for update day one, whichDay = 1.
    %  Note: whichDay = 0, -1, etc. can be used to plot days before the first
    %  update, novel, or novel2 day).
    % sessionToPlot specifices which sessions to include for lap plots. 1 =
    %  laps from all sessions on a day (e.g., plot all laps from an original session
    %  and an update session), 1 = laps from original sessions only, 2 =
    %  laps from update sessions only, 3 = laps from novel sessions only, and 
    %  4 = laps from novel2 sessions only.
    
    for an = length(params.animals)
        plotByDateOrType = 2; %0 = by date, 1 = by original day number, 2 = by update day number, 3 = by novel day number, 4 = by novel2 day number
        plotDate = [2501103]; %define date to plot for each mouse (rows = mice, columns = days)
        whichDay = [1]; %1 = first original/update/novel/novel2 day, 0 = day before update/novel/novel2 day
        sessionToPlot = 2; %0 = plot all, 1 = only plot original sessions on a given day, 2 = only plot update session on a given day, 3 = only plot novel session on a given day, 4 = only plot novel2 session on a given day
        doSessionPlots = 0; %plot data across all sessions for this mouse
        doLapPlots = 1;  %plot data across all laps on a given day for this mouse
        plotLapByBlock = 1; %0 = plot all laps individually, 1 = plot by blocks of 5, 2 = plot by block of specific trials
        numTrPerBlock = params.numTrPerBlock;
        whichBlocks = 'all'; %used with plotLapByBlock = 1; options = 'all', 'firstLast'
        %determine indices for this animal
        anindex = allindex(allindex(:,1) == params.animals(an), :);%indices for this animal
        anindex = anindex(anindex(:,6) == 2,:);%indices for active sessions

        if plotByDateOrType == 0 %plotting by date using plotDate, whichDay irrelevant
            thisp = size(plotDate,2);
            for p = 1:thisp
                plotBehaviorLapAndSession_DC(anindex,an,dirs,params,plotByDateOrType,plotDate(an,p),whichDay,sessionToPlot,doSessionPlots,doLapPlots,plotLapByBlock,numTrPerBlock,whichBlocks)
            end
        elseif plotByDateOrType == 1 || plotByDateOrType == 2 || plotByDateOrType == 3 || plotByDateOrType == 4 %plotting by type using whichDay, plotDate irrelevant
            thisp = length(whichDay);
            for p = 1:thisp
                plotBehaviorLapAndSession_DC(anindex,an,dirs,params,plotByDateOrType,plotDate(1),whichDay(p),sessionToPlot,doSessionPlots,doLapPlots,plotLapByBlock,numTrPerBlock,whichBlocks)
            end
        end%if plotByDateOrType == 0
    end%an

    %%%%%%%%%%%%%%%%
    %%%%% ROCs %%%%%
    %%%%%%%%%%%%%%%%

    for id = 1:length(params.rocID)%1 = speed, 2 = lickrate, 3 = deltalickrate
        %plot the AUC results
        dirROC = dir([dirs.saveoutputstructs, 'Data\Behavior\ROC']);
        ROCtoPlotInd = find(contains({dirROC.name}, ['ROC' '_' params.rocID{id}]),1,'last');%most recent data
        ROCtoPlot = load([dirROC(1).folder '\' dirROC(ROCtoPlotInd(1)).name]);
        ROC = ROCtoPlot.ROC;
        doROCPlots = 0; %ROC plots across sessions/animals, separated by track
        doAUCIndPlots = 1; %AUC plots across sessions for each animal
        doAUCGroupUpdatePlots = 0;%AUC plots across original sessions and update sessions by control vs. experimental group
        plotBehaviorROC_DC(allindex, uniqSess, dirs, ROC, params.rocID{id}, params, doROCPlots, doAUCIndPlots, doAUCGroupUpdatePlots)
    end

end%if plotBehavior

%% Gather neural data %%
if gatherNeuralData

    for i = 1:size(allindex,1) %loop through every session%need to remove some of my loops within this lol
        if allindex(i,5) > 0 %recording sessions only

            %%%%% session info %%%%%
            subj = [params.iden num2str(allindex(i,1))];
            sessDate = num2str(allindex(i,2));
            sessNum = num2str(allindex(i,3));
            trackInfo = num2str(allindex(i,4));
            %virmenDataPath = fullfile(dirs.virmenrawdata, ['\' subj '_' sessDate '_' sessNum '\' 'dataWithLickometer.mat']);
            virmenSessDataPath = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData\' subj '\' sessDate '_' sessNum '_' trackInfo]);
            neuralRawDataPath = fullfile(dirs.rawdata, [subj '_' sessDate]);
            processedDataPath = fullfile(dirs.processeddata,[subj '_' sessDate '\']);
            if ~exist(processedDataPath, 'dir')
                mkdir(processedDataPath);
            end
            saveNeuralPath = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' subj '\' sessDate '_' sessNum '_' trackInfo]);
            postSpikeSort = 1;

            if isfolder(virmenSessDataPath)%ensure behavior structs created using code above

                %%%%% make save path if not already created %%%%%
                if ~isfolder(saveNeuralPath)
                    mkdir(saveNeuralPath)
                end

                %%%%% create single ap bin file for all recordings on a given day %%%%%
                if ~isfolder([processedDataPath '\kilosort4'])
                    sprintf('Running Kilosort4 for %s_%s', subj, sessDate)
                    if params.iden == 'DC'
                        files = getfilenums(neuralRawDataPath);
                        getKilosort4Out_intan(subj, sessDate, neuralRawDataPath, dirs, files, params)%for Intan
                    else
                        getKilosort4Out(subj, sessDate, neuralRawDataPath, dirs)%for Neuropixels
                    end
                end

                if postSpikeSort

                    %%%%% create clusters, clustermetrics, and clusters_allrec structs %%%%%
                    if ~isfile([processedDataPath '\kilosort4\clusters.mat']) || params.rewrite.clusters
                        sprintf('Creating clusters structure for %s_%s_%s', subj, sessDate, sessNum)
                        Kilosort4_Pipeline_intan_DC(subj, sessDate, processedDataPath)%for Neuropixels
                    end

                    %%%%% create cell_metrics, noiseLevel, session, and spikes structs %%%%%
                    if ~isfile([processedDataPath '\kilosort4\' sprintf('%s_%s.cell_metrics.cellinfo.mat', subj, sessDate)]) || params.rewrite.cell_metrics
                        sprintf('Running CellExplorer for %s_%s_%s', subj, sessDate, sessNum)
                        addpath(genpath('\\ad.gatech.edu\bme\labs\singer\Danielle\code\CellExplorer'));
                        CellExplorer_Pipeline_DC(subj, sessDate, processedDataPath, params)%for Neuropixels
                    end

                    %%%%% add neural to behavior structs %%%%%
                    %Note: Gets session, lap, and trial stucts. Currently only
                    % includes complete laps and their trials.
                    if ~isfile([saveNeuralPath '\' 'rawDataBySessionNeural.mat']) || params.rewrite.neuralStructs
                        sprintf('Getting neural structs for %s_%s_%s', subj, sessDate, sessNum)
                        getNeuralStructs_linearDC(subj, sessDate, sessNum, params, virmenSessDataPath, processedDataPath, saveNeuralPath)
                    end

                    %%%%% get cell yield info for this session %%%%%
                    %Note: Uses rawDataBySessionNeural and only do session 3
                    if ~isfile([saveNeuralPath '\' 'cellYield.mat']) || params.rewrite.cellYield
                            sprintf('Getting cell yield info for %s_%s_%s', subj, sessDate, sessNum)
                            getCellYieldInfo(saveNeuralPath)
                    end

                    %%%%% get pyramidal layer info for this session %%%%%
                    %Note: Uses rawDataBySessionNeural and clusters_allrec structs
                    if ~isfile([saveNeuralPath '\' 'sessionPyrLayerInfo.mat']) || params.rewrite.pyrLayer
                        sprintf('Getting pyramidal layer info for %s_%s_%s', subj, sessDate, sessNum)
                        plotPyrLayer = 1;
                        selectManually = 1;
                        getPyrLayerInfo(subj, sessDate, sessNum, dirs, params, neuralRawDataPath, saveNeuralPath, plotPyrLayer, selectManually)
                    end

                    %%%%% get ripples for this session %%%%%
                    %Note: Uses rawDataBySessionNeural struct and chooses channel with most ripples that pass criteria
                    if ~isfield([saveNeuralPath '\' 'rawDataBySessionNeural'], 'ripplesGood') || params.rewrite.ripples
                        sprintf('Getting ripples for %s_%s_%s', subj, sessDate, sessNum)
                        plotRipples = 1;
                        getRipplesTmp_DC(dirs, params, saveNeuralPath, plotRipples)
                    end


                end%if postSpikeSort

            end%if isfile(virmenSessDataPath)

        end%if allindex(i,5) > 0
    end%i

end%if gatherNeuralData

%% Decoding %%
if doDecoding

%     for i = 1:size(allindex,1) %loop through every session
%         if allindex(i,5) > 0 %recording sessions only
% 
% 
% 
%             %%%%% session info %%%%%
%             subj = [params.iden num2str(allindex(i,1))];
%             sessDate = num2str(allindex(i,2));
%             sessNum = num2str(allindex(i,3));
%             trackInfo = num2str(allindex(i,4));
%             virmenSessDataPath = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData\' subj '\' sessDate '_' sessNum '_' trackInfo]);
%             neuralRawDataPath = fullfile(dirs.rawdata, [subj '_' sessDate]);
%             saveNeuralPath = fullfile(dirs.saveoutputstructs, ['Data\Neural\sessionData\' subj '\' sessDate '_' sessNum '_' trackInfo]);
% 
% 
% 
% 
end







