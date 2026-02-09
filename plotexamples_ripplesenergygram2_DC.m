function [mnP, maxP, allmntbdratio] = plotexamples_ripplesenergygram2_DC(savefigsdir, ripple, ripples, eeg, timearoundrip, freqnumerator,...
    freqdenominator, freqinterest, interactive, subj)
%based on plotexamples_ripplesenergygram2 2/3/26
%inputs needed
%   savefigsdir - i do think i want this automated
%   timearoundrip - to know how much to plot
%   freqnum and freqdenom - from params anyway
%   freqinterest - i assume
%   plot esamples - idk
%   interactive - so i can see as we go
%   eeg, ripple, ripples: remove thetas, tdbratio
%   we do not want to cycle through files, just one file
%updated ALP 3/18/21 to include brainRegion for title in plots

figure(3)
ratiolabel = [num2str(freqnumerator(1)), '-', num2str(freqnumerator(2)), '/', num2str(freqdenominator(1)), '-', num2str(freqdenominator(2)) ] ;
for w = 1:2
    mnP{w} = [];
    maxP{w} = [];
    subplot(3,1,w)
    hold on
    xlabel(['mean Power ratio: ', ratiolabel])
    ylabel(['max Power ratio: ', ratiolabel])
end
subplot(3,1,3)
hold on
xlabel(['mean Power ratio: ', ratiolabel])
ylabel('mean tdbratio')
subplot(3,1,1)
title('100msec around ripple')
subplot(3,1,2)
title([num2str(timearoundrip), ' sec around ripple'])
subplot(3,1,3)
title('100msec around ripple for ration, 2 sec around ripple for  tdbratio')
subtitle([num2str(subj), ' ', 'CA1'])
    
    %get sampling rates
    samprate = round(ripple.samprate);
    
    %make times vector
    a = max(size(ripple.data))/samprate; %ALP 4/1/21 updated to be max, though I think it was working the other way
    times = [0:1/samprate:a]';
    times = times(1:end-1);
    
    %get ripple envelope
    renv = ripple.data(:,3);
    smoothing_width = 0.004;
    kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
    renv = smoothvect(renv, kernel);
    
    %get ripple indices + time around
    rips = ripples.midind;
    ripperiods = [rips-timearoundrip*samprate rips+timearoundrip*samprate]; %take time before and after ripplem so ~2 seconds todal
    
    while any(ripperiods) & ripperiods(1,1)<0 %if first ripple too close to start, exclude it
        ripperiods = ripperiods(2:end,:);
    end
    while any(ripperiods) & ripperiods(end,2) > size(times,1) %if last ripple too close to end, exclude it
        ripperiods = ripperiods(1:end-1,:);
    end
    
    if any(ripperiods)
        outripperiods = [];
        for r = 1:size(ripperiods,1)
            %get equivalent length periods outside of rip periods
            endindex = ripperiods(r,2);
            if r == size(ripperiods,1)
                maxindex =  size(ripple.data);
            else
                maxindex = ripperiods(r+1,1);
            end
            while endindex+(samprate*timearoundrip*2) < maxindex %while long enough period between end of last period and start of next ripple period
                if all(renv(endindex:endindex+(samprate*timearoundrip*2)) < 3*ripples.std+ripples.baseline)%no SWR > 2std
                    outripperiods = [outripperiods; endindex endindex+(samprate*timearoundrip*2)];
                end
                endindex = endindex+(samprate*timearoundrip*2)+1;
                
            end
        end
        periods{2} = outripperiods ;
        periods{1} = ripperiods;
        
        %some plotting params
        limstd = ripples.baseline+ ripples.std*7;
        periodtype{1} = 'rip'; periodtype{2} = 'nonrip';
        
        %get number of per to plot - show about 20 (?) alp added 3/18/21
        if size(periods{1},1) > 20
            plotinterval = floor(size(periods{1},1)/20);
        else
            plotinterval = 1;
        end
        
        
        %plot LFP and intracell for each SWR
        for p =1%:2
            meantdb = nan(size(periods{p},1), 1);
            for r = 1:plotinterval:size(periods{p},1)
                plotind = [periods{p}(r,1): periods{p}(r,2)]; %indices to select for plotting, includes time around ripplss and combine overlapping periods
                psdind = rips(r)-0.050*samprate: rips(r)+0.050*samprate; %100msec centered on ripple mid index
                
                wcplotind = plotind*10-9;
                if max(plotind)<size(times,1)
                    %calculate energygram
                    % [energy phase freq energygramtime] = calcenergygram(lfp_filtered(plotind),  samprate, freqinterest(1), freqinterest(2), 10, 2);
                    %[energy phase freq time] = calcenergygram(data,  srate, lowf, highf, bandwidth, freqstep)
                    
                    %compute psd
                    for w = 1:2 %for short or long windows
                        if w == 1
                            [Pxx{w},F{w}] = pwelch(detrend(eeg.data(psdind)),[],[],[],samprate);
                        elseif w == 2
                            [Pxx{w},F{w}] = pwelch(detrend(eeg.data(plotind)),[],[],[],samprate);
                        end
                        freqs{w} = find(F{w}>=freqinterest(1) & F{w}<=freqinterest(2));
                        
                        %compute max & mean power in 150-250 and 250-350
                        subfreqs{1} = find(F{w}>=freqnumerator(1) & F{w}<=freqnumerator(2)); %ripple band
                        subfreqs{2} = find(F{w}>freqdenominator(1) & F{w}<=freqdenominator(2)); %above ripple, prob movement
                        for f = 1:2 %for 150-250 or 250-400
                            tempmnP{w}(1,f) = mean(Pxx{w}(subfreqs{f}));
                            tempmaxP{w}(1,f) = max(Pxx{w}(subfreqs{f}));
                        end
                        mnP{w} = [mnP{w}; tempmnP{w}]; %[numerator denominator]
                        maxP{w} = [maxP{w}; tempmaxP{w}];
                    end
                    
                    %mean tdbrateion
                    % meantdb(r) = mean(tdbratio.data(plotind));
                    % 
                    % if  plotexamples == 1 %(rem(r, 20) == 0 | tempmnP{1}(1)/tempmnP{1}(2) < 7) &&  plotexamples == 1 %plot every X ripple and ripples < 7 mean ratio
                    %     %plot
                    %     figure(1)
                    %     subplot(4,1,1)
                    %     plot(times(plotind), eeg.data(plotind), 'k') %eeg
                    %     xlim([times(plotind(1)) times(plotind(end))])
                    %     title('LFP')
                    % 
                    %     subplot(4,1,2)
                    %     plot(times(plotind), ripple.data((plotind),1), 'k') %rippleband
                    %     hold on
                    %     plot(times(plotind), renv((plotind),1), 'g') %smoothed envelope
                    %     plot([times(plotind(1)) times(plotind(end))], ripples.threshold*[1 1], '--')
                    %     plot([times(plotind(1)) times(plotind(end))], ripples.threshold*[-1 -1], '--')
                    %     ylim([-limstd limstd])
                    %     xlim([times(plotind(1)) times(plotind(end))])
                    %     title([ 'ripples: filtered 150-250Hz and envelope(green).  stdev(dash)',num2str(ripples{end}{end}{end}.std)])
                    % 
                    %     subplot(4,1,3)%power spectral density (mean)
                    %     for w = 1:2
                    %         plot(F{w}(freqs{w}), Pxx{w}(freqs{w}))
                    %     end
                    %     xlim(freqinterest)
                    %     %title(['psd. mean P ratio:', num2str(tempmnP{w}(1)/tempmnP{w}(2)), ', max ratio P ', num2str(tempmaxP{w}(1)/tempmaxP(2))])
                    % 
                    %     %plot tdbratio - remove this, do speed instead.
                    %     %will need to convert lfp indices to vr indices.
                    %     %need to figure that out
                    %     subplot(4,1,4)
                    %     hold on
                    %     plot(times(plotind), tdbratio.data((plotind),1), 'g', 'linewidth', 3) %theta/deltaratio
                    %     plot([times(plotind(1)) times(plotind(end))], thetas.threshold*[1 1], ':')
                    %     plot([times(plotind(1)) times(plotind(end))], (thetas.baseline-thetas.std)*[1 1], '--')
                    %     plot([times(plotind(1)) times(plotind(end))], thetas.baseline*[1 1], 'k-')
                    %     xlim([times(plotind(1)) times(plotind(end))])
                    %     title(['theta delta ratio (green), baseline (solid), 1std below baseline (dashed), peak threshold (dotted)'])
                    %     %subtitle(num2str(index));
                    %     subtitle(['index: ', num2str(index(1)),' ', num2str(index(2)), ' ', num2str(index(3)), periodtype{p},  ' period ', num2str(r), ' ', brainReg]);
                    % 
                    %     %                     figure(2)
                    %     %                     imagesc(energygramtime, freq, energy )
                    %     %                     colorbar('SouthOutside' )
                    %     %                     subtitle(['index: ', num2str(index), periodtype{p},  ' period ', num2str(r)])
                    %     %
                        if interactive
                            pause
                        else
                            ripplefigdir = fullfile(savefigsdir, 'RipplePeriods', filesep);
                            if ~exist(ripplefigdir); mkdir(ripplefigdir); end;
                            figfilename = [ripplefigdir 'rippleexamples' num2str(r)];
                            saveas(gcf, figfilename, 'png');
                        end
                    %end
                        %%%% i think this is the end of actual ripple
                        %%%% plotting stuff
                    %     figure(3)
                    %     for w = 1:2
                    %         subplot(3,1,w)
                    %         plot(tempmnP{w}(1)/tempmnP{w}(2), tempmaxP{w}(1)/tempmaxP{w}(2), '*r')
                    %     end
                    %     subplot(3,1,3)
                    %     plot(tempmnP{1}(1)/tempmnP{1}(2), meantdb(r), '*r')
                    % 
                    % 
                    %     clf(1)
                    %     %clf(2)
                    % end
                    % 
                    % if plotexamples == 1
                    %     figure(3)
                    %     for w = 1:2
                    %         subplot(3,1,w)
                    %         plot(tempmnP{w}(1)/tempmnP{w}(2), tempmaxP{w}(1)/tempmaxP{w}(2), '*b')
                    %     end
                    %     subplot(3,1,3)
                    %     plot(tempmnP{1}(1)/tempmnP{1}(2), meantdb(r), '*b')
                    %     subtitle(num2str(index));
                    % 
                    %     if interactive
                    %         pause
                    %     else
                    %         ripplefigdir = fullfile(savefigsdir, 'RipplePeriods', filesep);
                    %         if ~exist(ripplefigdir); mkdir(ripplefigdir); end;
                    %         figfilename = [ripplefigdir 'rippleenergygram' num2str(r)];
                    %         saveas(gcf, figfilename, 'png');
                    %     end
                    % end
                end
            end
            %numinclrips(i) = length(inclrips);
            %numexclrips(i) = length(exclrips);
            if ~isempty(periods{p})
                propripperiods(i) = sum(periods{p}(:,2)-periods{p}(:,1))/size(times,1);
            end
        end
        allmntbdratio = [allmntbdratio; meantdb(r)];
    end
    clear tdbratio ripples ripple eeg
    
end