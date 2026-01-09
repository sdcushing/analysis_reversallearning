function [statsByRestSession] = getRestBehaviorStats_linearJLK(rawDataBySession, virmen_fileInfo, saveBehaviorPath)
%adapted from getRestSessionStats_linearJLK

%% Get stats for rest session %%

%%%%% Calculate total time immobile %%%%%
immobility = find(rawDataBySession.rotVelo+rawDataBySession.transVelo==0);

immobilityEndInds = find(diff(immobility)>1);%last zero in string of no movement before movement

if ~isempty(immobilityEndInds)%if no movement detected all session, all values will be zeros
    
    immobilityLengths = [immobilityEndInds(1); diff(immobilityEndInds)];
    
    immobilityDur = nan(length(immobilityEndInds), 1);
    for ind = 1:length(immobilityEndInds)
        tmpimmobilityTimes = rawDataBySession.vrTime(immobility(immobilityEndInds(ind) - immobilityLengths(ind)+1:immobilityEndInds(ind)),1);
        tmpimmobilityDur = (tmpimmobilityTimes(end) - tmpimmobilityTimes(1));
        immobilityDur(ind) = tmpimmobilityDur;
    end
    
    statsByRestSession = [];
    statsByRestSession.sessionInfo = virmen_fileInfo.sessioninfo;
    statsByRestSession.totTimeImmobile = sum(immobilityDur);%seconds
    statsByRestSession.perTimeImmobile = (sum(immobilityDur) / ((rawDataBySession.vrTime(end) -  rawDataBySession.vrTime(1)))) * 100;%percent
    statsByRestSession.immobilityDurations = immobilityDur;%seconds
    statsByRestSession.immobilityMaxDuration = max(immobilityDur);%seconds
else
    statsByRestSession = [];
end


%% Save stats %%
filename = [saveBehaviorPath '\' 'statsByRestSession.mat'];
save(filename, 'statsByRestSession');

end%function