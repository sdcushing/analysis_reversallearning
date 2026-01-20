function [Azones, Rzones, NRzones, NevRzones] = getZoneInfo_linearJLK(fileInfo, subj)

%determine zones
if contains(subj, 'DC')
    Rzones = fileInfo.degReward;
    if strcmp(fileInfo.trackname, "TrackA'") || strcmp(fileInfo.sessioninfo, "DC21_251118_5")%update track
        NRzones = [36 206 326];
        NevRzones = [6 186 246];
            if sum(strcmp(subj, {'DC21'})) ~= 0
                NevRzones = [76,176,226];
            end
    elseif strcmp(fileInfo.trackname, "TrackA")%og track
        NRzones = [106 176 296];
        NevRzones = [NaN,NaN,NaN];
    end
elseif contains(subj, 'JK') || contains(subj, 'C')
    if strcmp(fileInfo.trackname, "TrackA'")%update track
        %non-reward zones = original anticipatory zones
        if sum(strcmp({'JK5' 'JK6'}, subj)) ~= 0
            NRzones = [40,130,230,310]; Rzones = [20,120,230,310];
        elseif sum(strcmp({'JK7' 'JK8' 'JK9' 'JK10'}, subj)) ~= 0
            NRzones = [40,130,230,310]; Rzones = [20,120,190,280];
        elseif sum(strcmp({'JK11' 'JK12'}, subj)) ~= 0
            NRzones = [40,130,230,310]; Rzones = [90,200,280];
        else
            NRzones = [20,130,230]; Rzones = [90,200,280];
        end
    
        %alternative non-reward zones = not original reward zones
        if sum(strcmp(subj, {'JK5', 'JK6'})) ~= 0
            NevRzones = [90,200,350];
        elseif sum(strcmp(subj, {'JK7', 'JK8', 'JK9', 'JK10'})) ~= 0
            NevRzones = [90,170,350];
        else
            NevRzones = [60,170,350];
        end
    else%original or novel track
        Rzones = fileInfo.degReward; NRzones = wrapTo360(Rzones + 3*fileInfo.cueSize); NevRzones = [NaN,NaN,NaN,NaN];
    end
elseif contains(subj, 'X')
    if isequal(fileInfo.degReward, [10, 150, 280])%update track
        %non-reward zones = original reward zones
        Rzones = fileInfo.degReward; % updated RZ
        NRzones = [50, 210, 330]-10; % original AZ
        %alternative non-reward zones = not original reward zones
        NevRzones = [80, 100, 240]; % NRzones from original track
    else%original or novel track
        Rzones = fileInfo.degReward; 
        NRzones = [80, 100, 240]; 
        NevRzones = [NaN,NaN,NaN];
    end
else
    sprintf('User must check getZoneInfo_linearJLK.m function.\n')
end

%anticipatory zones are always one zone before the reward zones
Azones = Rzones - fileInfo.cueSize;

end%function



