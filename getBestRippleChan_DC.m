function out = getBestRippleChan_DC(ripple, ripples)
%% this function finds the channel with the highest ripple power to be included for any data analyses
% based on getBestRippleChan_simple
% DC 01.30.26
% INPUTS:
%   ripple - ripple filtered lfp
%   ripples - actual detected ripples
% OUTPUTS:
%       top channel in terms of ripple power output from 0:31
% updated for 2 shank takahashi probes - ALP 11/22/19
% updated ALP 4/1/21 for a simple format for the preprocessing stream
    for i = 1:size(ripples, 2)%not sure what structure this variable is, will need to check
            chanIdx = i; 
            %get average ripple envelope
            rippleEnv = ripple.env(i,:);
            avgRipplePowertemp = mean(rippleEnv);
            
            %get peak ripple values
            if size(ripples(i).peak,1) > 0
                ripplePeaktemp = max(ripples(i).peak);
                rippleEnergytemp = max(ripples(i).energy);
            else
                ripplePeaktemp = nan;
                rippleEnergytemp = nan;
            end
        avgRipplePower(i,:) = avgRipplePowertemp'; %this value is derived from overall amplitude in the ripple band
        maxRipplePeak(i,:) = ripplePeaktemp; %this value is derived from actual ripple events
        maxRippleEnergy(i,:) = rippleEnergytemp; %this value is derived from actual ripple events
    end
    
    
    %         remove rows with missing files
    avgRipplePower = avgRipplePower(any(avgRipplePower,2),:);
    maxRipplePeak = maxRipplePeak(any(maxRipplePeak,2),:);
    maxRippleEnergy = maxRippleEnergy(any(maxRippleEnergy,2),:);
    
    %% plot these values as a sanity check
    %this section commented out for the simple version ALP 4/1/21
    %make your own version and enter the probe geometry  if you'd like
    %to plot these things!
    
    %         for chanIdx = 1:32
    %             avgRipplePower_geomorderedchan(1:size(avgRipplePower,1),chanIdx) = avgRipplePower(:,channelgeom.(probeside{p})(chanIdx)+1);
    %             maxRipplePeak_geomorderedchan(1:size(maxRipplePeak,1),chanIdx) = maxRipplePeak(:,channelgeom.(probeside{p})(chanIdx)+1);
    %             maxRippleEnergy_geomorderedchan(1:size(maxRippleEnergy,1),chanIdx) = maxRippleEnergy(:,channelgeom.(probeside{p})(chanIdx)+1);
    %         end
    
    %         figure(1); hold on; subplot(1,2,p); hold on;
    %         plot(avgRipplePower_geomorderedchan','LineWidth',2)
    %         title(ports{p})
    %         xlabel('Channels from Top to Bottom','FontSize',20)
    %         xticks(1:32); xlim([0 33])
    %         ylabel('Mean SWR Envelope')
    %
    %
    %         figure(2); hold on; subplot(1,2,p); hold on;
    %         plot(maxRipplePeak_geomorderedchan','LineWidth',2)
    %         title(ports{p})
    %         xlabel('Channels from Top to Bottom','FontSize',20)
    %         xticks(1:32); xlim([0 33])
    %         ylabel('Max SWR Peak')
    %
    %
    %         figure(3); hold on; subplot(1,2,p); hold on
    %         plot(maxRippleEnergy_geomorderedchan','LineWidth',2)
    %         title(ports{p})
    %         xlabel('Channels from Top to Bottom','FontSize',20)
    %         xticks(1:32); xlim([0 33])
    %         ylabel('Max SWR Energy')
    %
    
    %% get ripple channel
    %highest average ripple power
    if size(avgRipplePower,1) > 1
        meanRipPower = nanmean(avgRipplePower);
    else
        meanRipPower = avgRipplePower;
    end
    [~, bestIdx] = max(avgRipplePower);
    bestRippleChan.channel = bestIdx;
    bestRippleChan.info = 'channel with largest mean power in the ripple band - eeg index'; 
    bestRippleChan.date = datestr(now); 
    st = dbstack; 
    bestRippleChan.callFun = {st(:).name}; 

out = bestRippleChan; 
end