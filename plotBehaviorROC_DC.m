function plotBehaviorROC_DC(allindex, uniqSess, dirs, ROC, rocID, params, doROCPlots, doAUCIndPlots, doAUCGroupUpdatePlots)

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Saving Info %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

figdir = fullfile([dirs.savefigures 'Behavior\ROC\']);
if ~isfolder(figdir)
    mkdir(figdir)
end
cd(figdir)
%save figure including info about which subjects were plotted
ids = []; for s = 1:length(ROC.sessInfo); thisid = ROC.sessInfo{s}(1); ids = [ids thisid]; end
ids = num2str(unique(ids));

if doROCPlots

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% ROC Plots by Track %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %original anticipatory vs control zones (original reward zones + 30 deg)
    figure
    hold on
    colors = [];
    colors = hsv(length(ROC.og_all.X));
    colorctr = 0;
    for r = 1:length(ROC.og_all.AUC)
        if sum(isnan(ROC.og_all.AUC))>=1%if any non-original track sessions
            if r < find(isnan(ROC.og_all.AUC),1,'first')
                plot(ROC.og_all.X{r}, ROC.og_all.Y{r},'Color',colors(r,:),'LineWidth', 2)
            elseif r > find(isnan(ROC.og_all.AUC),1,'first') && ~isnan(ROC.og_all.AUC(r))
                colorctr = colorctr + 1;
                plot(ROC.og_all.X{r}, ROC.og_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            end
        else%only original track sessions
            colorctr = colorctr + 1;
            plot(ROC.og_all.X{r}, ROC.og_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
        end
    end
    legend('Box','off')
    plot(0:.1:1,0:.1:1,'-k')
    title(sprintf('Original AZ vs Non RZ %s ROC: %s %s',rocID,params.iden,ids))
    figname = fullfile(figdir, sprintf('og_%s_ROC_across_sessions',rocID));
    print(gcf,figname,'-dpng','-r300')

    %update anticipatory vs original reward zones
    if isfield(ROC.up_all,'X')
        figure
        hold on
        colors = [];
        colors = hsv(length(ROC.up_all.X));
        colorctr = 0;
        newanctr = 0;
        for r = 1:length(ROC.up_all.AUC)
            if r >= find(~isnan(ROC.up_all.AUC),1,'first') && r <= length(ROC.up_all.AUC)/2
                colorctr = colorctr + 1;
                plot(ROC.up_all.X{r}, ROC.up_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            elseif r > length(ROC.up_all.AUC)/2 && ~isnan(ROC.up_all.AUC(r))
                newanctr = newanctr + 1;
                if newanctr == 1
                    colorctr = 0;
                end
                colorctr = colorctr + 1;
                plot(ROC.up_all.X{r}, ROC.up_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            end
        end
        legend('Box','off')
        plot(0:.1:1,0:.1:1,'-k')
        title(sprintf('Update AZ vs Original RZ %s ROC: %s %s',rocID,params.iden,ids))
        figname = fullfile(figdir, sprintf('up_%s_ROC_across_sessions',rocID));
        print(gcf,figname,'-dpng','-r300')
    end%if isfield(ROC.up_all,'X')

    %update anticipatory vs never rewarded zones
    if isfield(ROC.up_UAZvNevRZ,'X')
        figure
        hold on
        colors = [];
        colors = hsv(length(ROC.up_UAZvNevRZ.X));
        colorctr = 0;
        newanctr = 0;
        for r = 1:length(ROC.up_UAZvNevRZ.X)
            if r >= find(~isnan(ROC.up_UAZvNevRZ.AUC),1,'first') && r <= length(ROC.up_UAZvNevRZ.AUC)/2
                colorctr = colorctr + 1;
                plot(ROC.up_UAZvNevRZ.X{r}, ROC.up_UAZvNevRZ.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            elseif r > length(ROC.up_UAZvNevRZ.AUC)/2 && ~isnan(ROC.up_UAZvNevRZ.AUC(r))
                newanctr = newanctr + 1;
                if newanctr == 1
                    colorctr = 0;
                end
                colorctr = colorctr + 1;
                plot(ROC.up_UAZvNevRZ.X{r}, ROC.up_UAZvNevRZ.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            end
        end
        legend('Box','off')
        plot(0:.1:1,0:.1:1,'-k')
        title(sprintf('Update AZ vs Never RZ %s ROC: %s %s',rocID,params.iden,ids))
        figname = fullfile(figdir, sprintf('up_UAZvNevRZ_%s_ROC_across_sessions',rocID));
        print(gcf,figname,'-dpng','-r300')
    end

    %original reward vs never rewarded zones
    if isfield(ROC.up_ORZvNevRZ,'X')
        figure
        hold on
        colors = [];
        colors = hsv(length(ROC.up_ORZvNevRZ.X));
        colorctr = 0;
        newanctr = 0;
        for r = 1:length(ROC.up_ORZvNevRZ.X)
            if r >= find(~isnan(ROC.up_ORZvNevRZ.AUC),1,'first') && r <= length(ROC.up_ORZvNevRZ.AUC)/2
                colorctr = colorctr + 1;
                plot(ROC.up_ORZvNevRZ.X{r}, ROC.up_ORZvNevRZ.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            elseif r > length(ROC.up_ORZvNevRZ.AUC)/2 && ~isnan(ROC.up_ORZvNevRZ.AUC(r))
                newanctr = newanctr + 1;
                if newanctr == 1
                    colorctr = 0;
                end
                colorctr = colorctr + 1;
                plot(ROC.up_ORZvNevRZ.X{r}, ROC.up_ORZvNevRZ.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            end
        end
        legend('Box','off')
        plot(0:.1:1,0:.1:1,'-k')
        title(sprintf('Update Original RZ Vs. Never Rewarded Zones %s ROC: %s %s',rocID,params.iden,ids))
        figname = fullfile(figdir, sprintf('up_ORZvNevRZ_%s_ROC_across_sessions',rocID));
        print(gcf,figname,'-dpng','-r300')
    end

    %novel vs control zones (novel reward zones + 30 deg)
    if isfield(ROC.nov_all,'X')
        figure
        hold on
        colors = [];
        colors = hsv(length(ROC.nov_all.X));
        colorctr = 0;
        newanctr = 0;
        for r = 1:length(ROC.nov_all.AUC)
            if r >= find(~isnan(ROC.nov_all.AUC),1,'first') && r <= length(ROC.nov_all.AUC)/2
                colorctr = colorctr + 1;
                plot(ROC.nov_all.X{r}, ROC.nov_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            elseif r > length(ROC.nov_all.AUC)/2 && ~isnan(ROC.nov_all.AUC(r))
                newanctr = newanctr + 1;
                if newanctr == 1
                    colorctr = 0;
                end
                colorctr = colorctr + 1;
                plot(ROC.nov_all.X{r}, ROC.nov_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            end
        end
        legend('Box','off')
        plot(0:.1:1,0:.1:1,'-k')
        title(sprintf('Novel Zones %s ROC: %s %s',rocID,params.iden,ids))
        figname = fullfile(figdir, sprintf('nov_%s_ROC_across_sessions',rocID));
        print(gcf,figname,'-dpng','-r300')
    end%if isfield(ROC.nov_all,'X')

    %novel2 vs control zones (novel2 reward zones + 30 deg)
    if isfield(ROC.nov2_all,'X')
        figure
        hold on
        colors = [];
        colors = hsv(length(ROC.nov2_all.X));
        colorctr = 0;
        newanctr = 0;
        for r = 1:length(ROC.nov2_all.AUC)
            if r >= find(~isnan(ROC.nov2_all.AUC),1,'first') && r <= length(ROC.nov2_all.AUC)/2
                colorctr = colorctr + 1;
                plot(ROC.nov2_all.X{r}, ROC.nov2_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            elseif r > length(ROC.nov2_all.AUC)/2 && ~isnan(ROC.nov2_all.AUC(r))
                newanctr = newanctr + 1;
                if newanctr == 1
                    colorctr = 0;
                end
                colorctr = colorctr + 1;
                plot(ROC.nov2_all.X{r}, ROC.nov2_all.Y{r},'Color',colors(colorctr,:),'LineWidth', 2)
            end
        end
        legend('Box','off')
        plot(0:.1:1,0:.1:1,'-k')
        title(sprintf('Novel2 Zones %s ROC: %s %s',rocID,params.iden,ids))
        figname = fullfile(figdir, sprintf('nov2_%s_ROC_across_sessions',rocID));
        print(gcf,figname,'-dpng','-r300')
    end%if isfield(ROC.nov2_all,'X')

end%if doROCPlots

if doAUCIndPlots

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% AUC Plot Across Tracks: Individual %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %star used for single day, line used for multiple days
    lgdNames = '';
    AUCfig = figure('Name','AUC');
    hold on

    %original anticipatory zones vs control zones (original reward zones + 30 deg)
    if length(ROC.og_all.AUC) == 1
        plot(ROC.og_all.AUC, '*r')
        lgdNames = [lgdNames {'OAZvNonRZ'}];
    elseif length(ROC.og_all.AUC) > 1
        plot(ROC.og_all.AUC, '-r')
        lgdNames = [lgdNames {'OAZvNonRZ'}];
    end

    %update anticipatory zones vs original reward zones
    if sum(~isnan(ROC.up_all.AUC)) == 1
        plot(ROC.up_all.AUC, '*g')
        lgdNames = [lgdNames {'UAZvORZ'}];
    elseif sum(~isnan(ROC.up_all.AUC)) > 1
        plot(ROC.up_all.AUC, '-g')
        lgdNames = [lgdNames {'UAZvORZ'}];
        %deal with separated single days
        tmpDiffs = find(diff(~isnan([ROC.up_all.AUC nan])));
        tmpDiffsToPlot = [];
        for d = 1:length(tmpDiffs)-1
            if tmpDiffs(d+1) - tmpDiffs(d) == -2
                tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
            end
        end
        if ~isempty(tmpDiffsToPlot)
            plot(tmpDiffsToPlot,ROC.up_all.AUC(tmpDiffsToPlot), '*g')
            lgdNames = [lgdNames {'UAZvORZ'}];
        end
    end

    %update anticipatory zones vs never rewarded control zones
    if sum(~isnan(ROC.up_UAZvNevRZ.AUC)) == 1
        plot(ROC.up_UAZvNevRZ.AUC, '*b')
        lgdNames = [lgdNames {'UAZvNevRZ'}];
    elseif sum(~isnan(ROC.up_UAZvNevRZ.AUC)) > 1
        plot(ROC.up_UAZvNevRZ.AUC, '-b')
        lgdNames = [lgdNames {'UAZvNevRZ'}];
        %deal with separated single days
        tmpDiffs = find(diff(~isnan([ROC.up_UAZvNevRZ.AUC nan])));
        tmpDiffsToPlot = [];
        for d = 1:length(tmpDiffs)-1
            if tmpDiffs(d+1) - tmpDiffs(d) == -2
                tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
            end
        end
        if ~isempty(tmpDiffsToPlot)
            plot(tmpDiffsToPlot,ROC.up_UAZvNevRZ.AUC(tmpDiffsToPlot), '*b')
            lgdNames = [lgdNames {'UAZvNevRZ'}];
        end
    end

    %original reward zones vs never rewarded control zones
    if sum(~isnan(ROC.up_ORZvNevRZ.AUC)) == 1
        plot(ROC.up_ORZvNevRZ.AUC, '*', 'Color', [1 .6 .6])
        lgdNames = [lgdNames {'ORZvNevRZ'}];
    elseif sum(~isnan(ROC.up_ORZvNevRZ.AUC)) > 1
        plot(ROC.up_ORZvNevRZ.AUC,  '-', 'Color', [1 .6 .6])
        lgdNames = [lgdNames {'ORZvNevRZ'}];
        %deal with separated single days
        tmpDiffs = find(diff(~isnan([ROC.up_ORZvNevRZ.AUC nan])));
        tmpDiffsToPlot = [];
        for d = 1:length(tmpDiffs)-1
            if tmpDiffs(d+1) - tmpDiffs(d) == -2
                tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
            end
        end
        if ~isempty(tmpDiffsToPlot)
            plot(tmpDiffsToPlot,ROC.up_ORZvNevRZ.AUC(tmpDiffsToPlot), '*', 'Color', [1 .6 .6])
            lgdNames = [lgdNames {'ORZvNevRZ'}];
        end
    end

    %novel anticipatory zones vs control zones (novel reward zones + 30 deg)
    if sum(~isnan(ROC.nov_all.AUC)) == 1
        plot(ROC.nov_all.AUC, '*m')
        lgdNames = [lgdNames {'NovAZvNonRZ'}];
    elseif sum(~isnan(ROC.nov_all.AUC)) > 1
        plot(ROC.nov_all.AUC, '-m')
        lgdNames = [lgdNames {'NovAZvNonRZ'}];
        %deal with separated single days
        tmpDiffs = find(diff(~isnan([ROC.nov_all.AUC nan])));
        tmpDiffsToPlot = [];
        for d = 1:length(tmpDiffs)-1
            if tmpDiffs(d+1) - tmpDiffs(d) == -2
                tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
            end
        end
        if ~isempty(tmpDiffsToPlot)
            plot(tmpDiffsToPlot,ROC.nov_all.AUC(tmpDiffsToPlot), '*m')
            lgdNames = [lgdNames {'NovAZvNonRZ'}];
        end
    end

    %novel2 anticipatory zones vs control zones (novel2 reward zones + 30 deg)
    if sum(~isnan(ROC.nov2_all.AUC)) == 1
        plot(ROC.nov2_all.AUC, '*y')
        lgdNames = [lgdNames {'Nov2AZvNonRZ'}];
    elseif sum(~isnan(ROC.nov2_all.AUC)) > 1
        plot(ROC.nov2_all.AUC, '-y')
        lgdNames = [lgdNames {'Nov2AZvNonRZ'}];
        %deal with separated single days
        tmpDiffs = find(diff(~isnan([ROC.nov2_all.AUC nan])));
        tmpDiffsToPlot = [];
        for d = 1:length(tmpDiffs)-1
            if tmpDiffs(d+1) - tmpDiffs(d) == -2
                tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
            end
        end
        if ~isempty(tmpDiffsToPlot)
            plot(tmpDiffsToPlot,ROC.nov2_all.AUC(tmpDiffsToPlot), '*y')
            lgdNames = [lgdNames {'Nov2AZvNonRZ'}];
        end
    end

    %add number label for each subject
    lbls = [];
    for r = 1:length(ROC.sessInfo)
        lbls = [lbls ROC.sessInfo{r}(1)];
    end
    a = unique(lbls);
    tmpxtick = zeros(length(a),1);
    for an = 1:length(a)
        %find the number of sessions for this animal, divide by two to find the
        %middle session, then add that value to the first session for the
        %animal
        tmpxtick(an) = ceil( (find(lbls==a(an),1,'last') - find(lbls==a(an),1,'first')) / 2 ) + find(lbls==a(an),1,'first');
    end
    xticks(tmpxtick)
    xticklabels(a)
    xlabel([params.iden '#'])
    ylabel('AUC')
    ylim([0 1])
    title(sprintf('%s AUC across Animals',rocID))

    %loop through sessions to mark 'x' on sessions below threshold
    for r = 1:length(ROC.sessInfo)
        dayOgTrNum = 0;
        dayUpTrNum = 0;
        dayNovTrNum = 0;
        dayNov2TrNum = 0;
        for s = 1:size(ROC.sessInfo{r},1)
            %load session info
            sessionInfo = ROC.sessInfo{r}(s,:);
            f1 = fullfile(dirs.saveoutputstructs, ['Data\Behavior\sessionData\' params.iden num2str(sessionInfo(1))], ...
                [num2str(sessionInfo(2)) '_' num2str(sessionInfo(3)) '_' num2str(sessionInfo(4))], 'statsByLap.mat');
            lapData = load(f1);
            lapData = lapData.statsByLap;

            %save number of trials separately for original and update tracks
            if ~isempty(lapData)
                if ROC.sessInfo{r}(s,4) == 1
                    dayOgTrNum = dayOgTrNum + length(lapData.numRewards);
                elseif ROC.sessInfo{r}(s,4) == 2
                    dayUpTrNum = dayUpTrNum + length(lapData.numRewards);
                elseif ROC.sessInfo{r}(s,4) == 3
                    dayNovTrNum = dayNovTrNum + length(lapData.numRewards);
                elseif ROC.sessInfo{r}(s,4) == 4
                    dayNov2TrNum = dayNov2TrNum + length(lapData.numRewards);
                end%session type
            end%if ~isempty(lapData)
        end%s

        %mark an x on days with fewer than 4 trials
        if dayOgTrNum < 4
            figure(AUCfig)
            hold on
            plot(r, ROC.og_all.AUC(r),'xk')
        end
        if dayUpTrNum < 4
            figure(AUCfig)
            hold on
            plot(r, ROC.up_all.AUC(r),'xk')
            plot(r, ROC.up_UAZvNevRZ.AUC(r),'xk')
            plot(r, ROC.up_ORZvNevRZ.AUC(r),'xk')
        end
        if dayNovTrNum < 4
            figure(AUCfig)
            hold on
            plot(r, ROC.nov_all.AUC(r),'xk')
        end
        if dayNov2TrNum < 4
            figure(AUCfig)
            hold on
            plot(r, ROC.nov2_all.AUC(r),'xk')
        end

    end%sessions

    %legend
    lgdNames = [lgdNames {'< 4 Trials'}];
    legend(lgdNames,'Box','off','Location','southeast','FontSize',6)

    %save
    figname = fullfile(figdir, sprintf('%s_AUC_across_sessions',rocID));
    print(gcf,figname,'-dpng','-r300')

end%if doAUCIndPlots

if doAUCGroupUpdatePlots

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% AUC Plot Across Tracks: Group Average %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %loop though groups, using only experimental mice that have been unblineded
    ogDataPerAn = nan(2,10,2);%group x animals x AUC vals
    upDataPerAn = nan(2,10,4);%group x animals x AUC vals
    upUAZvNevRZDataPerAn = nan(2,10,4);%group x animals x AUC vals
    upORZvNevRZDataPerAn = nan(2,10,4);%group x animals x AUC vals
    gNames = {'APOE3', 'APOE4'};
    for g = 1:2%1 = control, 2 = experimental
        if g == 1%control (APOE3)
            mice2Use = ismember(params.animals,params.controlGroup) & ismember(params.animals,params.testedMice);
        elseif g == 2%experimental (APOE4)
            mice2Use = ismember(params.animals,params.experimentalGroup)& ismember(params.animals,params.testedMice);
        end%if g == 1
        upDay = nan(1,sum(mice2Use));
        upFinalSess = nan(1,sum(mice2Use));
        if sum(mice2Use) > 0
            miceNums = params.animals(mice2Use);
            for an = 1:length(miceNums)
                %loop through sessions to find update session
                for r = 1:length(ROC.sessInfo)
                    for s = 1:size(ROC.sessInfo{r},1)
                        sessionInfo = ROC.sessInfo{r}(s,:);
                        if sessionInfo(1) == miceNums(an) && sessionInfo(4) == 2 && isnan(upDay(an))%update session for mouse in this group
                            upDay(an) = r;
                        end%if
                        if sessionInfo(1) == miceNums(an) && sessionInfo(4) == 2
                            upFinalSess(an) = r;
                        end%if
                    end%s
                end%r

                %save data per animal for this group
                ogDataPerAn(g,an,:) = ROC.og_all.AUC(upDay(an)-1:upDay(an));
                upDataPerAn(g,an,1:upFinalSess(an)-upDay(an)+1) = ROC.up_all.AUC(upDay(an):upFinalSess(an));
                upUAZvNevRZDataPerAn(g,an,1:upFinalSess(an)-upDay(an)+1) = ROC.up_UAZvNevRZ.AUC(upDay(an):upFinalSess(an));
                upORZvNevRZDataPerAn(g,an,1:upFinalSess(an)-upDay(an)+1) = ROC.up_ORZvNevRZ.AUC(upDay(an):upFinalSess(an));
            end%an
        end%if sum(mice2use) > 0
    end%g
                  
    %find average and SEM across animals within each group; output = group x day
    ogDataAvg = squeeze(nanmean(ogDataPerAn,2));
    ogDataSEM = squeeze(nanstd(ogDataPerAn,1,2)) ./ squeeze(sqrt(sum(~isnan(ogDataPerAn),2)));
    upDataAvg = squeeze(nanmean(upDataPerAn,2));
    upDataSEM = squeeze(nanstd(upDataPerAn,1,2)) ./ squeeze(sqrt(sum(~isnan(upDataPerAn),2)));
    upUAZvNevRZDataAvg = squeeze(nanmean(upUAZvNevRZDataPerAn,2));
    upUAZvNevRZDataSEM = squeeze(nanstd(upUAZvNevRZDataPerAn,1,2)) ./ squeeze(sqrt(sum(~isnan(upUAZvNevRZDataPerAn),2)));
    upORZvNevRZDataAvg = squeeze(nanmean(upORZvNevRZDataPerAn,2));
    upORZvNevRZDataSEM = squeeze(nanstd(upORZvNevRZDataPerAn,1,2)) ./ squeeze(sqrt(sum(~isnan(upORZvNevRZDataPerAn),2)));

    %plot separately by session type
    for dayType = 1:2%1 = original, 2 = update
        figure; hold on
        lgdNames = [];
        plotDataAvg = [];
        plotDataSEM = [];
        p = [];
        if dayType == 1%original
            colors = [1 .6 .6; 1 0 0];%contorl = light red, experimental = red
            plotDataAvg = ogDataAvg;
            plotDataSEM = ogDataSEM;
            %plot by group
            for g = 1:2%1 = control, 2 = experimental
                p(g) = plot(plotDataAvg(g,:)', '-', 'Color', colors(g,:), 'LineWidth',1.5);
                lgdNames = [lgdNames {[gNames{g} ' ' 'OAZvNonRZ']}];
                errorbar([plotDataAvg(g,:)' plotDataAvg(g,:)'], [plotDataSEM(g,:)' plotDataSEM(g,:)'], 'Color', colors(g,:), 'LineWidth',1)
            end
            xlim([0.5 2.5])
            xticks([1 2])
            xticklabels([1 2])
            xlabel('Day')
            ylabel('AUC')
            ylim([0.3 1])
            legend(p, lgdNames,'Box','off','Location','southeast','FontSize',6)
            title(sprintf('%s AUC by Group: Original Sessions',rocID))
            figname = fullfile(figdir, sprintf('%s_AUC_by_group_original_sessions',rocID));
        elseif dayType == 2%update
            plotCtr = 0;
            vNames = {'UAZvORZ','UAZvNevRZ','ORZvNevRZ'};
            for vars = 1:3%1 = update AZ vs. og RZ, 2 = update AZ vs never rewarded zone, 3 = og RZ vs. never rewarded zone
                if vars == 1%update AZ vs. og RZ
                    colors = [.6 1 .6; 0 1 0];%contorl = light green, experimental = green
                    plotDataAvg = upDataAvg;
                    plotDataSEM = upDataSEM;
                elseif vars == 2%update AZ vs never rewarded zone
                    colors = [.6 .6 1; 0 0 1];%contorl = light blue, experimental = blue
                    plotDataAvg = upUAZvNevRZDataAvg;
                    plotDataSEM = upUAZvNevRZDataSEM;
                elseif vars == 3%og RZ vs. never rewarded zone
                    colors = [1 .6 .6; 1 0 0];%contorl = light red, experimental = red
                    plotDataAvg = upORZvNevRZDataAvg;
                    plotDataSEM = upORZvNevRZDataSEM;
                end
                %plot by group
                for g = 1:2%1 = control, 2 = experimental
                    plotCtr = plotCtr + 1;
                    p(plotCtr) = plot(plotDataAvg(g,:)', '-', 'Color', colors(g,:), 'LineWidth',1.5);
                    lgdNames = [lgdNames {[gNames{g} ' ' vNames{vars}]}];
                    errorbar([plotDataAvg(g,:)' plotDataAvg(g,:)'], [plotDataSEM(g,:)' plotDataSEM(g,:)'], 'Color', colors(g,:), 'LineWidth',1)
                end%g
            end%vars
            xlim([0.5 4.5])
            xticks([1 2 3 4])
            xticklabels([1 2 3 4])
            xlabel('Day')
            ylabel('AUC')
            ylim([0 1])
            legend(p, lgdNames,'Box','off','Location','southeast','FontSize',6)
            title(sprintf('%s AUC by Group: Update Sessions',rocID))
            figname = fullfile(figdir, sprintf('%s_AUC_by_group_update_sessions',rocID));
        end%if dayType == 1

        %save
        print(gcf,figname,'-dpng','-r300')
    end%dayType

end%if doAUCGroupUpdatePlots


% % % %
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%%%% AUC Plot Update Session First5Last5 %%%%%
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %
% % % % %%%%% plot first and last 5 trials of the update day %%%%%
% % % % %loop through animals to find update day 2 session index
% % % % upDayInd = [];
% % % % for an = 1:length(unique(uniqSess(:,1)))
% % % %     allAn = unique(uniqSess(:,1));
% % % %     thisAn = allAn(an);
% % % %     tmpUpSess = find(allindex(:,1) == thisAn & allindex(:,4) == 2 & allindex(:,6) == 2,3,'first');
% % % %     upDayInd(an) = find(uniqSess(:,1) == allindex(tmpUpSess(2),1) & uniqSess(:,2) == allindex(tmpUpSess(2),2));
% % % % end
% % % %
% % % % %anticipatory zones vs original reward zones
% % % % %block1 and 2 as a line
% % % % figure('Name','AUC');
% % % % hold on
% % % % xTickNames = '';
% % % % for an = 1:length(unique(uniqSess(:,1)))
% % % %     plot((an-1:an)+an , [ROC.up_block1.AUC(upDayInd(an)) ROC.up_block2.AUC(upDayInd(an))], '-ok','MarkerFaceColor','k')
% % % %     xTickNames{an} = sprintf('JK%d',allAn(an));
% % % % end
% % % % xlim([0.5 length(allAn)*2+.5])
% % % % xticks([1.5:2:length(allAn)*2])
% % % % xticklabels(xTickNames)
% % % % ylabel('AUC')
% % % % ylim([.3 1])
% % % % title(sprintf('%s AUC: First 5 and Last 5 Trials\n of Update Day 2 across Animals',rocID))
% % % %
% % % % %save
% % % % figname = fullfile(figdir, sprintf('%s_AUC_across_sessions_first5last5',rocID));
% % % % print(gcf,figname,'-dpdf','-r300')

% % % % CURRENTLY COMMENTING OUT UAZVSNEVRZ AND OAZVSNEVRZ
% % % % %star used for single day, line used for multiple days
% % % % lgdNames = '';
% % % % AUCfig = figure('Name','AUC');
% % % % hold on
% % % % %update anticipatory zones vs never rewarded control zones
% % % % %block1
% % % % color = "#80B3FF";
% % % % if sum(~isnan(ROC.upnev_block1.AUC)) == 1
% % % %     plot(ROC.upnev_block1.AUC, '*','Color',color)
% % % %     lgdNames = [lgdNames {'UpAlt_First5'}];
% % % % elseif sum(~isnan(ROC.upnev_block1.AUC)) > 1
% % % %     plot(ROC.upnev_block1.AUC, '-','Color',color)
% % % %     lgdNames = [lgdNames {'UpAlt_First5'}];
% % % %     %deal with separated single days
% % % %     tmpDiffs = find(diff(~isnan([ROC.upnev_block1.AUC nan])));
% % % %     tmpDiffsToPlot = [];
% % % %     for d = 1:length(tmpDiffs)-1
% % % %         if tmpDiffs(d+1) - tmpDiffs(d) == 1
% % % %             tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
% % % %         end
% % % %     end
% % % %     if ~isempty(tmpDiffsToPlot)
% % % %         plot(tmpDiffsToPlot,ROC.upnev_block1.AUC(tmpDiffsToPlot), '*','Color',color)
% % % %         lgdNames = [lgdNames {'UpAlt_First5'}];
% % % %     end
% % % % end
% % % % %block2
% % % % if sum(~isnan(ROC.upnev_block2.AUC)) == 1
% % % %     plot(ROC.upnev_block2.AUC, '*b')
% % % %     lgdNames = [lgdNames {'UpAlt_Last5'}];
% % % % elseif sum(~isnan(ROC.upnev_block2.AUC)) > 1
% % % %     plot(ROC.upnev_block2.AUC, '-b')
% % % %     lgdNames = [lgdNames {'UpAlt_Last5'}];
% % % %     %deal with separated single days
% % % %     tmpDiffs = find(diff(~isnan([ROC.upnev_block2.AUC nan])));
% % % %     tmpDiffsToPlot = [];
% % % %     for d = 1:length(tmpDiffs)-1
% % % %         if tmpDiffs(d+1) - tmpDiffs(d) == 1
% % % %             tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
% % % %         end
% % % %     end
% % % %     if ~isempty(tmpDiffsToPlot)
% % % %         plot(tmpDiffsToPlot,ROC.upnev_block2.AUC(tmpDiffsToPlot), '*b')
% % % %         lgdNames = [lgdNames {'UpAlt_Last5'}];
% % % %     end
% % % % end
% % % %
% % % % %original reward zones vs never rewarded control zones
% % % % %block1
% % % % color = [1 .7 .7];
% % % % if sum(~isnan(ROC.upnev2_block1.AUC)) == 1
% % % %     plot(ROC.upnev2_block1.AUC, '*','Color',color)
% % % %     lgdNames = [lgdNames {'UpAlt2_First5'}];
% % % % elseif sum(~isnan(ROC.upnev2_block1.AUC)) > 1
% % % %     plot(ROC.upnev2_block1.AUC, '-','Color',color)
% % % %     lgdNames = [lgdNames {'UpAlt2_First5'}];
% % % %     %deal with separated single days
% % % %     tmpDiffs = find(diff(~isnan([ROC.upnev2_block1.AUC nan])));
% % % %     tmpDiffsToPlot = [];
% % % %     for d = 1:length(tmpDiffs)-1
% % % %         if tmpDiffs(d+1) - tmpDiffs(d) == 1
% % % %             tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
% % % %         end
% % % %     end
% % % %     if ~isempty(tmpDiffsToPlot)
% % % %         plot(tmpDiffsToPlot,ROC.upnev2_block1.AUC(tmpDiffsToPlot), '*','Color',color)
% % % %         lgdNames = [lgdNames {'UpAlt2_First5'}];
% % % %     end
% % % % end
% % % % %block2
% % % % if sum(~isnan(ROC.upnev2_block2.AUC)) == 1
% % % %     plot(ROC.upnev2_block2.AUC, '*r')
% % % %     lgdNames = [lgdNames {'UpAlt2_Last5'}];
% % % % elseif sum(~isnan(ROC.upnev2_block2.AUC)) > 1
% % % %     plot(ROC.upnev2_block2.AUC, '-r')
% % % %     lgdNames = [lgdNames {'UpAlt2_Last5'}];
% % % %     %deal with separated single days
% % % %     tmpDiffs = find(diff(~isnan([ROC.upnev2_block2.AUC nan])));
% % % %     tmpDiffsToPlot = [];
% % % %     for d = 1:length(tmpDiffs)-1
% % % %         if tmpDiffs(d+1) - tmpDiffs(d) == 1
% % % %             tmpDiffsToPlot = [tmpDiffsToPlot tmpDiffs(d+1)];
% % % %         end
% % % %     end
% % % %     if ~isempty(tmpDiffsToPlot)
% % % %         plot(tmpDiffsToPlot,ROC.upnev2_block2.AUC(tmpDiffsToPlot), '*r')
% % % %         lgdNames = [lgdNames {'UpAlt2_Last5'}];
% % % %     end
% % % % end
% % % %
% % % % %add number label for each subject
% % % % lbls = [];
% % % % for r = find(~isnan(ROC.upnev_block1.AUC),1,'first'):find(~isnan(ROC.upnev_block1.AUC),1,'last')
% % % %     lbls = [lbls ROC.sessInfo{r}(1)];
% % % % end
% % % % a = unique(lbls);
% % % % tmpxtick = zeros(length(a),1);
% % % % for an = 1:length(a)
% % % %     %find the number of sessions for this animal, divide by two to find the
% % % %     %middle session, then add that value to the first session for the
% % % %     %animal
% % % %     tmpxtick(an) = ceil( (find(lbls==a(an),1,'last') - find(lbls==a(an),1,'first')) / 2 ) + find(lbls==a(an),1,'first') + find(~isnan(ROC.upnev_block1.AUC),1,'first');
% % % % end
% % % % xticks(tmpxtick)
% % % % xticklabels(a)
% % % % xlabel([params.iden '#'])
% % % % ylabel('AUC')
% % % % ylim([.3 1])
% % % % title(sprintf('%s AUC First5Last5 across Animals',rocID))
% % % % %legend
% % % % legend(lgdNames,'Box','off','Location','southeast','FontSize',6, 'Interpreter', 'none')
% % % %
% % % % %save
% % % % figname = fullfile(figdir, sprintf('%s_AUC_across_sessions_first5last5',rocID));
% % % % print(gcf,figname,'-dpdf','-r300')

cd('\\ad.gatech.edu\bme\labs\singer\Danielle\code\AnalysisCode\Neuropixels_analyses');
end%function