function plotBehaviorLapAndSession_JLK(anindex,an,dirs,params,plotByDateOrType,plotDate,whichDay,sessionToPlot,doSessionPlots,doLapPlots,plotLapByBlock,numTrPerBlock, whichBlocks)
%adapted from script_AnnularTrackAnalysisJLK.m

%% plot behavioral data over sessions %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% initialize variables %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% determine first update + novel day and plot day %%%%%

%days with active session for this mouse
animalUniqSess = unique(anindex(:,2));
%first original, update, and novel day
ogDaySessionDate = anindex(find(anindex(:,4) == 1, 1, 'first'),2);
upDaySessionDate = anindex(find(anindex(:,4) == 2, 1, 'first'),2);
if ~isempty(upDaySessionDate)
    upDayInd = find(unique(anindex(:,2))==upDaySessionDate, 1, 'first');
end
novDaySessionDate = anindex(find(anindex(:,4) == 3, 1, 'first'),2);
if ~isempty(novDaySessionDate)
    novDayInd = find(unique(anindex(:,2))==novDaySessionDate, 1, 'first');
end
nov2DaySessionDate = anindex(find(anindex(:,4) == 4, 1, 'first'),2);
if ~isempty(nov2DaySessionDate)
    nov2DayInd = find(unique(anindex(:,2))==nov2DaySessionDate, 1, 'first');
end

%day to plot
plotDay = [];
if plotByDateOrType==0%by date
    plotDay = plotDate;
    plotDayInd = find(animalUniqSess==plotDate);
elseif plotByDateOrType==1%by original day number
    plotDay = animalUniqSess(find(animalUniqSess==ogDaySessionDate)+whichDay-1);
    plotDayInd = find(animalUniqSess==ogDaySessionDate)+whichDay-1;
elseif plotByDateOrType==2%by update day number
    plotDay = animalUniqSess(find(animalUniqSess==upDaySessionDate)+whichDay-1);
    plotDayInd = find(animalUniqSess==upDaySessionDate)+whichDay-1;
elseif plotByDateOrType==3%by novel day number
    plotDay = animalUniqSess(find(animalUniqSess==novDaySessionDate)+whichDay-1);
    plotDayInd = find(animalUniqSess==novDaySessionDate)+whichDay-1;
elseif plotByDateOrType==4%by novel2 day number
    plotDay = animalUniqSess(find(animalUniqSess==nov2DaySessionDate)+whichDay-1);
    plotDayInd = find(animalUniqSess==nov2DaySessionDate)+whichDay-1;
end

%%%%% initialize zones + variables to plot %%%%%
OGrZones = NaN(1,4);
UPrZones = NaN(1,4);
NOVrZones = NaN(1,4);
NOV2rZones = NaN(1,4);
avgLickRate = zeros(size(anindex,1), 360/params.binsize_deg);
avgSpeed = zeros(size(anindex,1), 360/params.binsize_deg);
tmptrLickRate = [];
trLickRate = [];
tmptrSpeed = [];
trSpeed = [];
dayType = [];
plotVsFirstDay = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% loop through sessions to collect data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:size(anindex,1)
    sessionInfo = anindex(i,:);
    f1 = fullfile(dirs.saveoutputstructs, '\Data\Behavior\sessionData\', [params.iden num2str(sessionInfo(1), '%01d')], ...
        [num2str(sessionInfo(2)) '_' num2str(sessionInfo(3)) '_' num2str(sessionInfo(4))], 'statsByLap.mat');
    f2 = fullfile(dirs.saveoutputstructs, '\Data\Behavior\sessionData\', [params.iden num2str(sessionInfo(1), '%01d')], ...
        [num2str(sessionInfo(2)) '_' num2str(sessionInfo(3)) '_' num2str(sessionInfo(4))], 'statsBySession.mat');

    if isfile(f1) && isfile(f2)
        lapData = load(f1);
        lapData = lapData.statsByLap;
        sessionData = load(f2);
        sessionData = sessionData.statsBySession;

        if ~isempty(lapData)
            %record first update and novel session and zones for plotting
            if strcmp(lapData.fileInfo.trackname, "TrackA")
                OGrZones = lapData.fileInfo.degReward;
            end
            if isempty(find(anindex(:,4) == 2, 1, 'first'))%no update day
                upDaySession = NaN;
                UPrZones = [NaN,NaN,NaN,NaN];
            elseif ~isempty(find(anindex(:,4) == 2, 1, 'first'))%update day
                upDaySession = find(anindex(:,4) == 2, 1, 'first');
                if strcmp(lapData.fileInfo.trackname, "TrackA'")
                    UPrZones = lapData.fileInfo.degReward;
                end
            end
            if isempty(find(anindex(:,4) == 3, 1, 'first'))%no novel day
                novDaySession = NaN;
                NOVrZones = [NaN,NaN,NaN,NaN];
            elseif ~isempty(find(anindex(:,4) == 3, 1, 'first'))%novel day
                novDaySession = find(anindex(:,4) == 3, 1, 'first');
                if strcmp(lapData.fileInfo.trackname, 'TrackB')
                    NOVrZones = lapData.fileInfo.degReward;
                end
            end
            if isempty(find(anindex(:,4) == 4, 1, 'first'))%no novel2 day
                nov2DaySession = NaN;
                NOV2rZones = [NaN,NaN,NaN,NaN];
            elseif ~isempty(find(anindex(:,4) == 4, 1, 'first'))%novel2 day
                nov2DaySession = find(anindex(:,4) == 4, 1, 'first');
                if strcmp(lapData.fileInfo.trackname, 'TrackC')
                    NOV2rZones = lapData.fileInfo.degReward;
                end
            end
    
            %ensure the session includes rewards
            if ~isempty(sessionData.lickBehavior)
                avgLickRate(i, :) = sessionData.lickBehavior.lickRateAvg;
                avgSpeed(i, :) = sessionData.velocBehavior.velocCountsSmoothSum;
    
                %collect some information from all laps for a given day
                if anindex(i,2) == plotDay
                    %gather data based on update/novel only or not
                    if sessionToPlot==1 && anindex(i,4)==1 %plot only original session for this day
                        sessionType = 'Original';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==2 && anindex(i,4)==2 %plot only update session for this day
                        sessionType = 'Update';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==3 && anindex(i,4)==3 %plot only novel session for this day
                        sessionType = 'Novel';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==4 && anindex(i,4)==4 %plot only novel2 session for this day
                        sessionType = 'Novel2';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==1 && sum(anindex(anindex(:,2)==plotDay,4)==1)==0 %accidentally only plotting original sessions when animal never completed one on this day;  auto plot all laps
                        sessionType = 'All';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==2 && sum(anindex(anindex(:,2)==plotDay,4)==2)==0 %accidentally only plotting update sessions when animal never completed one on this day;  auto plot all laps
                        sessionType = 'All';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==3 && sum(anindex(anindex(:,2)==plotDay,4)==3)==0 %accidentally only plotting novel sessions when animal never completed one on this day; auto plot all laps
                        sessionType = 'All';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==4 && sum(anindex(anindex(:,2)==plotDay,4)==4)==0 %accidentally only plotting novel2 sessions when animal never completed one on this day; auto plot all laps
                        sessionType = 'All';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    elseif sessionToPlot==0 %plot all sessions for this day
                        sessionType = 'All';
                        tmptrLickRate = lapData.lickRateSmooth;
                        trLickRate = [trLickRate; tmptrLickRate];
                        tmptrSpeed = lapData.velocCountsSmooth;
                        trSpeed = [trSpeed; tmptrSpeed];
                    end
    
                    %record session type
                    if anindex(i,4) == 1
                        dayType = 'Original';
                        plotVsFirstDay = plotDayInd; %day to plot vs. first original day (e.g., original day 1)
                    elseif anindex(i,4) == 2
                        dayType = 'Update';
                        plotVsFirstDay = plotDayInd - upDayInd + 1; %day to plot vs. first update day (e.g., update day 1)
                    elseif anindex(i,4) == 3
                        dayType = 'Novel';
                        plotVsFirstDay = plotDayInd - novDayInd + 1; %day to plot vs. first novel day (e.g., novel day 1)
                    elseif anindex(i,4) == 4
                        dayType = 'Novel2';
                        plotVsFirstDay = plotDayInd - nov2DayInd + 1; %day to plot vs. first novel day (e.g., novel day 1)
                    end
    
                end%if anindex(i,2) == plotDay
            end%if ~isempty(sessionData.lickBehavior)
        end%if ~isempty(lapData)
    end%if isfile(f1) && isfile(f2)
end%i

if doSessionPlots

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Session Plots %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %Avg Lick Rate Across Sessions, Linegraph
    figure; clear ax
    colors = cbrewer('div','RdBu', size(anindex,1)+3);
    colors = flip(colors);
    ax = axes('NextPlot','add','Box','off','XLim',[0 360/params.binsize_deg]);
    img = arrayfun( @(x) plot(avgLickRate(x,:), 'Color', abs(colors(x+3,:)), 'DisplayName', num2str(x)), 1:size(anindex,1));
    vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g--');
    vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r--');
    vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm--');
    vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y--');
    xticks(0: 40 / params.binsize_deg: 360 / params.binsize_deg)
    xticklabels(0:40:360)
    xlabel(ax, 'Position (bin)')
    ylabel(ax, 'Average lick rate (licks/s)')
    title(ax, [params.iden num2str(params.animals(an)) ' Average Lick Rate Over Sessions'])
    legend show
    legend boxoff
    figdir = fullfile([dirs.savefigures 'Behavior\Example_TrainingDays\' params.iden num2str(params.animals(an))]);
    if ~isfolder(figdir)
        mkdir(figdir)
    end
    cd(figdir)
    figname = fullfile(figdir, 'lick_rate_over_position_across_sessions');
    print(gcf,figname,'-dpng','-r300')

    %Avg Lick Rate Across Sessions, Heatmap
    figure; clear ax
    colormap(flipud(gray));
    imagesc(avgLickRate); box off;
    c = colorbar; c.Label.String = 'Lick rate (licks/s)';
    vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g--');
    vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r--');
    vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm--');
    vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y--');
    hline(upDaySession, 'g--')
    hline(novDaySession, 'm--')
    hline(nov2DaySession, 'y--')
    xticks(0: 40 / params.binsize_deg: 360 / params.binsize_deg)
    xticklabels(0:40:360)
    xlabel('Position (bin)'); ylabel('Session');
    title([params.iden num2str(params.animals(an)) ' Average Lick Rate Over Sessions'])
    figdir = fullfile([dirs.savefigures 'Behavior\Example_TrainingDays\' params.iden num2str(params.animals(an))]);
    if ~isfolder(figdir)
        mkdir(figdir)
    end
    figname = fullfile(figdir, 'lick_rate_heatmap');
    print(gcf,figname,'-dpng','-r300')

    %Avg Speed Across Sessions, Linegraph
    figure; clear ax
    colors = cbrewer('div','RdBu', size(anindex,1)+3);
    colors = flip(colors);
    ax = axes('NextPlot','add','Box','off','XLim',[0 360/params.binsize_deg]);
    img = arrayfun( @(x) plot(avgSpeed(x,:), 'Color', abs(colors(x+3,:)), 'DisplayName', num2str(x)), 1:size(anindex,1));
    vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g--');
    vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r--');
    vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm--');
    vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y--');
    xticks(0: 40 / params.binsize_deg: 360 / params.binsize_deg)
    xticklabels(0:40:360)
    xlabel(ax, 'Position (bin)')
    ylabel(ax, 'Average speed (degs/s)')
    title(ax, [params.iden num2str(params.animals(an)) ' Average Speed Over Sessions'])
    legend show
    legend boxoff
    figname = fullfile(figdir, 'mean_speed_over_position_across_sessions');
    print(gcf,figname,'-dpng','-r300')

    %Avg Lick Rate Across Sessions, Heatmap
    figure; clear ax
    colormap(flipud(gray));
    imagesc(avgSpeed); box off;
    c = colorbar; c.Label.String = 'Speed (deg/s)';
    vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g--');
    vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r--');
    vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm--');
    vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y--');
    hline(upDaySession, 'g--')
    hline(novDaySession, 'm--')
    hline(nov2DaySession, 'y--')
    xticks(0: 40 / params.binsize_deg: 360 / params.binsize_deg)
    xticklabels(0:40:360)
    xlabel('Position (bin)'); ylabel('Session');
    title([params.iden num2str(params.animals(an)) ' Average Speed Over Sessions'])
    figdir = fullfile([dirs.savefigures 'Behavior\Example_TrainingDays\' params.iden num2str(params.animals(an))]);
    if ~isfolder(figdir)
        mkdir(figdir)
    end
    figname = fullfile(figdir, 'mean_speed_heatmap');
    print(gcf,figname,'-dpng','-r300')
end%if doSessionPlots

if doLapPlots

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Plot Laps of a Given Day %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Lick Rate Across Laps, Linegraph
    figure
    if plotLapByBlock == 1%plot blocks of laps
        hold on
        numBlocks = ceil(size(trLickRate,1)/numTrPerBlock);
        if strcmp(whichBlocks, 'all')
            useBlocks = 1:numBlocks;
        elseif strcmp(whichBlocks, 'firstLast')
            useBlocks = [1, numBlocks];
        end%if strcmp(whichBlocks, 'all')
        colors = cbrewer('qual','Set2', length(useBlocks)+3);
        lgdNames = '';
        for bl = useBlocks
            trLickRateTmp = (bl-1)*numTrPerBlock+1:(bl-1)*numTrPerBlock+1+numTrPerBlock-1;
            trLickRateTmp = trLickRateTmp(trLickRateTmp <= size(trLickRate,1));%ensure last block only goes to last trial number
            tmpMn = mean(trLickRate(trLickRateTmp,:),1);
            plot(tmpMn, 'Color',colors(bl,:))
            lgdNames = [lgdNames {sprintf('Block%d Mn',bl)}];
            if length(trLickRateTmp) > 1
                tmpSEM = std(trLickRate(trLickRateTmp,:),1) / sqrt(length(trLickRateTmp));
                fill([1:size(trLickRate,2), fliplr(1:size(trLickRate,2))], [tmpMn+tmpSEM, fliplr(tmpMn-tmpSEM)], colors(bl,:) ,'FaceAlpha', .18, 'EdgeColor', 'none');
                lgdNames = [lgdNames {sprintf('Block%d SEM',bl)}];
            end%if length(trLickRateTmp) > 1
        end%bl
        legend(lgdNames,'Box','off','Location','northeast','FontSize',6,'NumColumns',2)
        titleName = sprintf('%s%d Lick Rate Over %d Lap Blocks for %s Day %d (%d)', params.iden, params.animals(an), numTrPerBlock, dayType, plotVsFirstDay, plotDay);
        figName = sprintf('lick_rate_over_position_across_lap_blocks_%d_%s', plotDay, sessionType);
    elseif plotLapByBlock == 2%plot blocks of laps of specific trials
        hold on
        useTrials = [1:numTrPerBlock; size(trLickRate,1)-numTrPerBlock+1:size(trLickRate,1)]; %for now, first n and last n trials, where n = num trials in a block
        colors = [0,0,1; 0,0,0];
        lgdNames = '';
        for bl = 1:size(useTrials,1)
            tmpMn = mean(trLickRate(useTrials(bl,:),:),1);
            plot(tmpMn, 'Color',colors(bl,:))
            lgdNames = [lgdNames {sprintf('Block%d Mn',bl)}];
            tmpSEM = std(trLickRate(useTrials(bl,:),:),1) / sqrt(size(useTrials,2));
            fill([1:size(trLickRate,2), fliplr(1:size(trLickRate,2))], [tmpMn+tmpSEM, fliplr(tmpMn-tmpSEM)], colors(bl,:) ,'FaceAlpha', .18, 'EdgeColor', 'none');
            lgdNames = [lgdNames {sprintf('Block%d SEM',bl)}];
        end%bl
        legend(lgdNames,'Box','off','Location','northeast','FontSize',6,'NumColumns',2)
        titleName = sprintf('%s%d Lick Rate Over %d Lap Blocks for %s Day %d (%d)', params.iden, params.animals(an), numTrPerBlock, dayType, plotVsFirstDay, plotDay);
        figName = sprintf('lick_rate_over_position_across_lap_blocks_spectr_%d_%s', plotDay, sessionType);
    elseif plotLapByBlock == 0%plot all laps individually
        hold on
        colors = cbrewer('qual','Set2', size(trLickRate,1)+3);
        img = arrayfun( @(x) plot(trLickRate(x,:), 'Color', abs(colors(x+3,:)), 'DisplayName', num2str(x)), 1:size(trLickRate,1));
        legend('Box','off','Location','northeast','FontSize',6,'NumColumns',2)
        titleName = sprintf('%s%d Lick Rate Over Laps for %s Day %d (%d)', params.iden, params.animals(an), dayType, plotVsFirstDay, plotDay);
        figName = sprintf('lick_rate_over_position_across_laps_%d_%s', plotDay, sessionType);
    end%plotByBlock
    ylim([0 7])
    xticks(0: 40 / params.binsize_deg: 360 / params.binsize_deg)
    xticklabels(0:40:360)
    xlabel('Position (bin)')
    ylabel('Lick rate (licks/s)')
    title(titleName)
    %plot zones
    if strcmp(sessionType, 'Original') || strcmp(sessionType, 'Update')
        vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g-');%update reward zone
        vline([UPrZones-10] ./ params.binsize_deg, 'g--');%update anticipatory zone
        vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zone
        vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
    elseif strcmp(sessionType, 'Novel')
        vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm-');%novel reward zone
        vline([NOVrZones-10] ./ params.binsize_deg, 'm--');%novel anticipatory zone
    elseif strcmp(sessionType, 'Novel2')
        vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y-');%novel2 reward zone
        vline([NOV2rZones-10] ./ params.binsize_deg, 'y--');%novel2 anticipatory zone
    elseif strcmp(sessionType, 'All') %for All sessions, draw zones based on day type
        if strcmp(dayType, 'Original')
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zone
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        elseif strcmp(dayType, 'Update')
            vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g-');%update reward zone
            vline([UPrZones-10] ./ params.binsize_deg, 'g--');%update anticipatory zone
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zone
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        elseif strcmp(dayType, 'Novel')
            vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm-');%novel reward zone
            vline([NOVrZones-10] ./ params.binsize_deg, 'm--');%novel anticipatory zone
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zon
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        elseif strcmp(dayType, 'Novel2')
            vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y-');
            vline([NOV2rZones-10] ./ params.binsize_deg, 'y--');%novel2 anticipatory zone
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%assuming original trials on a novel2 day
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        end
    end%plot zones
    figdir = fullfile([dirs.savefigures 'Behavior\Example_TrainingDays\' params.iden num2str(params.animals(an))]);
    if ~isfolder(figdir)
        mkdir(figdir)
    end
    figname = fullfile(figdir, figName);
    print(gcf,figname,'-dpng','-r300')

    %Speed Across Laps, Linegraph
    figure
    if plotLapByBlock == 1%plot blocks of laps
        hold on
        numBlocks = ceil(size(trSpeed,1)/numTrPerBlock);
        if strcmp(whichBlocks, 'all')
            useBlocks = 1:numBlocks;
        elseif strcmp(whichBlocks, 'firstLast')
            useBlocks = [1, numBlocks];
        end%if strcmp(whichBlocks, 'all')
        colors = cbrewer('qual','Set2', length(useBlocks)+3);
        lgdNames = '';
        for bl = useBlocks
            trLickRateTmp = (bl-1)*numTrPerBlock+1:(bl-1)*numTrPerBlock+1+numTrPerBlock-1;
            trLickRateTmp = trLickRateTmp(trLickRateTmp <= size(trSpeed,1));%ensure last block only goes to last trial number
            tmpMn = mean(trSpeed(trLickRateTmp,:),1);
            plot(tmpMn, 'Color',colors(bl,:))
            lgdNames = [lgdNames {sprintf('Block%d Mn',bl)}];
            if length(trLickRateTmp) > 1
                tmpSEM = std(trSpeed(trLickRateTmp,:),1) / sqrt(length(trLickRateTmp));
                fill([1:size(trSpeed,2), fliplr(1:size(trSpeed,2))], [tmpMn+tmpSEM, fliplr(tmpMn-tmpSEM)], colors(bl,:) ,'FaceAlpha', .18, 'EdgeColor', 'none');
                lgdNames = [lgdNames {sprintf('Block%d SEM',bl)}];
            end%if length(trLickRateTmp) > 1
        end%bl
        legend(lgdNames,'Box','off','Location','northeast','FontSize',6,'NumColumns',2)
        titleName = sprintf('%s%d Speed Over %d Lap Blocks for %s Day %d (%d)', params.iden, params.animals(an), numTrPerBlock, dayType, plotVsFirstDay, plotDay);
        figName = sprintf('mean_speed_over_position_across_lap_blocks_%d_%s', plotDay, sessionType);
    elseif plotLapByBlock == 2%plot blocks of laps of specific trials
        hold on
        useTrials = [1:numTrPerBlock; size(trSpeed,1)-numTrPerBlock+1:size(trSpeed,1)]; %for now, first n and last n trials, where n = num trials in a block
        colors = [0,0,1; 0,0,0];
        lgdNames = '';
        for bl = 1:size(useTrials,1)
            tmpMn = mean(trSpeed(useTrials(bl,:),:),1);
            plot(tmpMn, 'Color',colors(bl,:))
            lgdNames = [lgdNames {sprintf('Block%d Mn',bl)}];
            tmpSEM = std(trSpeed(useTrials(bl,:),:),1) / sqrt(size(useTrials,2));
            fill([1:size(trSpeed,2), fliplr(1:size(trSpeed,2))], [tmpMn+tmpSEM, fliplr(tmpMn-tmpSEM)], colors(bl,:) ,'FaceAlpha', .18, 'EdgeColor', 'none');
            lgdNames = [lgdNames {sprintf('Block%d SEM',bl)}];
        end%bl
        legend(lgdNames,'Box','off','Location','northeast','FontSize',6,'NumColumns',2)
        titleName = sprintf('%s%d Speed Over %d Lap Blocks for %s Day %d (%d)', params.iden, params.animals(an), numTrPerBlock, dayType, plotVsFirstDay, plotDay);
        figName = sprintf('mean_speed_over_position_across_lap_blocks_spectr_%d_%s', plotDay, sessionType);
    elseif plotLapByBlock == 0%plot all laps individually
        hold on
        colors = cbrewer('qual','Set2', size(trSpeed,1)+3);
        img = arrayfun( @(x) plot(trSpeed(x,:), 'Color', abs(colors(x+3,:)), 'DisplayName', num2str(x)), 1:size(trSpeed,1));
        legend('Box','off','Location','northeast','FontSize',6,'NumColumns',2)
        titleName = sprintf('%s%d Speed Over Laps for %s Day %d (%d)', params.iden, params.animals(an), dayType, plotVsFirstDay, plotDay);
        figName = sprintf('mean_speed_over_position_across_laps_%d_%s', plotDay, sessionType);
    end%plotByBlock
    ylim([0 floor(max(trSpeed,[],'all'))+2])
    xticks(0: 40 / params.binsize_deg: 360 / params.binsize_deg)
    xticklabels(0:40:360)
    xlabel('Position (bin)')
    ylabel('Speed (degs/s)')
    title(titleName)
    legend show
    legend boxoff
    %plot zones
    if strcmp(sessionType, 'Original') || strcmp(sessionType, 'Update')
        vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g-');%update reward zone
        vline([UPrZones-10] ./ params.binsize_deg, 'g--');%update anticipatory zone
        vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zone
        vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
    elseif strcmp(sessionType, 'Novel')
        vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm-');%novel reward zone
        vline([NOVrZones-10] ./ params.binsize_deg, 'm--');%novel anticipatory zone
    elseif strcmp(sessionType, 'Novel2')
        vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y-');%novel2 reward zone
        vline([NOV2rZones-10] ./ params.binsize_deg, 'y--');%novel2 anticipatory zone
    elseif strcmp(sessionType, 'All') %for All sessions, draw zones based on day type
        if strcmp(dayType, 'Original')
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zone
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        elseif strcmp(dayType, 'Update')
            vline([UPrZones, UPrZones + 10] ./ params.binsize_deg, 'g-');%update reward zone
            vline([UPrZones-10] ./ params.binsize_deg, 'g--');%update anticipatory zone
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zone
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        elseif strcmp(dayType, 'Novel')
            vline([NOVrZones, NOVrZones + 10] ./ params.binsize_deg, 'm-');%novel reward zone
            vline([NOVrZones-10] ./ params.binsize_deg, 'm--');%novel anticipatory zone
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%original reward zon
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        elseif strcmp(dayType, 'Novel2')
            vline([NOV2rZones, NOV2rZones + 10] ./ params.binsize_deg, 'y-');
            vline([NOV2rZones-10] ./ params.binsize_deg, 'y--');%novel2 anticipatory zone
            vline([OGrZones, OGrZones + 10] ./ params.binsize_deg, 'r-');%assuming original trials on a novel2 day
            vline([OGrZones-10] ./ params.binsize_deg, 'r--');%original anticipatory zone
        end
    end%plot zones
    figname = fullfile(figdir, figName);
    print(gcf,figname,'-dpng','-r300')
end%if doLapPlots


end%function
