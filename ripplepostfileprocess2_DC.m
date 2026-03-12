function ripplepostfileprocess2_DC(ripples, tdbratio, eeg, timearoundrip, freqnumerator, freqdenominator, ...
    ratiothresh, outlierindices_allchan, currentDeg, vrTime, lfpTime, varargin)
% based on ripplepostfileprocess2, compatible with new pipeline
%inputs
%   ripples - actual ripples (mid ind., samprate, startind, endind, )
%   lfp data - filtered (eeg- data, dounsample)
%   thetas - theta periods (samprate, baseline)
%   tdbratio - ratio theta delta beta (samprate, data)
%   'exclude' , 1 or 0, default 0
%           1: exclude outlier ripples
% INPUTS
%   Ratiothresh -  peaks over peakthresh will be considered outlier ripples
% OUTPUTS
%   excluded-      list of excluded ripples per epoch: [animal date cell file #outlier ripples.thresh]
%
% run after running extractripples to exclude events if > ratio threash  mean power
% ratio (150-250/250-400) for 400msec around rip to exclude events that are
% likely due to movement noise (>200Hz) instead of SWR

% assign the options
exclude = 0;
applyspeed = [];
applyMUA = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'exclude'
            exclude = varargin{option+1};
        case 'applyspeed'
            applyspeed = varargin{option+1};
        case 'applyMUA'
            applyMUA = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

excluded = [];
% go through each file in flist and filter it
    for fnum = 1:size(ripples, 2)%loop through channels
        tdbcurration = tdbratio.data((fnum),:);
        if any(ripples(fnum).midind)
            
            %get ratio of freq power
            %get this function, update with no index DC
            [ratio] = findpowerratioripplevsabove2(ripples(fnum), eeg(fnum,:), timearoundrip, freqnumerator, freqdenominator);
            if length(ratio) ~= length(ripples(fnum).midind)
                error('ratio is not being calculated properly in ripple post processing')
            end
            excl1 = ratio<ratiothresh; %ripples that should b excluded
            incl1 = ratio>=ratiothresh; %ripples that should be included
            
            %find if meet tdbratio criteria as well
            rips = ripples(fnum).midind;
            ripperiods = [rips-0.5*ripples(fnum).samprate rips+0.5*ripples(fnum).samprate]; %take time before and after ripplem so 1 seconds todal
            ripperiods(ripperiods(:,1)<1,1) = 1; %if start before start of recording, set to 1
            ripperiods(ripperiods(:,2)>length(tdbcurration), 2) = length(tdbcurration); %if start after end of recording set to end of recording
            
            for r = 1:size(ripperiods)
                meantdb(r,:) = mean(tdbratio.data(ripperiods(r,1):ripperiods(r,2)));%average tdbratio during 1 sec around center of rip
            end
            %             excl2 = meantdb>thetas{index(1)}{index(2)}{index(3)}.baselinethreshold; %exclude if average thetadeltabeta ratio is > thetas.baseline threshold
            %             incl2 = meantdb<=thetas{index(1)}{index(2)}{index(3)}.baselinethreshold;
            %%%% CHANGED TO BASELINE NOT BASELINE THRESH - 6.19.18 SP
            %removed for time being. think no longer necessary and movement
            %will be enough. DC 2.26.26
            %excl2 = meantdb>thetas.baseline; %exclude if average thetadeltabeta ratio is > thetas.baseline 
            %incl2 = meantdb<=thetas.baseline;
            
            %excl if outlier falls into ripple period
            %             tempout = getbinindex(eeg{end}{end}{end}.outlierindices, [ripples{index(1)}{index(2)}{index(3)}.startind ripples{index(1)}{index(2)}{index(3)}.endind]);
            tempout = getbinindex(outlierindices_allchan{fnum}, [ripples(fnum).startind ripples(fnum).endind]);
            excl3 = zeros(size(rips,1),1);
            excl3(unique(tempout(tempout~=0))) = 1;
            incl3 = ~excl3;
            
            %excl if max value is greater than maxvalthresh (default 1500)
            maxvalthresh = 1500;
            tempout2 = zeros(size(rips,1),1);
            for i = 1:length(ripples(fnum).startind)
                eegvalues = eeg(fnum,(ripples(fnum).startind(i):ripples(fnum).endind(i)));
                if any(eegvalues > maxvalthresh) || any(eegvalues < -maxvalthresh)
                    tempout2(i) = 1;
                end
            end
            excl4 = tempout2;
            incl4 = ~excl4;
            
            if isempty(excl1)
                excl1 = false(size(rips,1),1);
                incl1 = ~excl1;
            end
            
            %excl if zero MUA spike
            if applyMUA == 1
                spikeidx = spikes.spikeindex;
                downsamp = eeg.downsample;
                rips = [ripples(fnum).startind, ripples(fnum).endind];
                rips = rips.*downsamp; %turn it into spike index at 30000 samprate
                
                excl5 = false(size(rips,1),1);
                for ii = 1:size(rips,1)
                    excl5(ii) = sum(isExcluded(spikeidx, rips(ii,:))) == 0;
                end
                incl5 = ~excl5;
            else
                excl5 = false(size(rips,1),1);
                incl5 = ~excl5;
            end
            
            %excl if during period above speed threshold
            % load rawpos file to apply speed treshold
            if ~isempty(applyspeed)
                speedTh = applyspeed(1); %in deg/s, set on 01.05.2021
                t_sec = applyspeed(2); %in seconds, time range before and after ripple mid index
                %downsamp = eeg.downsample;%removing bc we have lfpTime
                %already made. changed downsamp to 1 in below equations
                rips = ripples(fnum).midind;
                %rips = rips.*downsamp; %turn it into spike index at 30000 samprate
                ripperiods = [rips-t_sec*ripples(fnum).samprate*1, rips+t_sec*ripples(fnum).samprate*1]; %take time before and after ripplem so 2 seconds total - had to make it 30k samp rate bc using spike index
                ripperiods(ripperiods(:,1)<1,1) = 1; %if start before start of recording, set to 1
                ripperiods(ripperiods(:,2)>length(tdbcurration)*1, 2) = length(tdbcurration)*1; %if start after end of recording set to end of recording
                
                excl6 = false(size(ripperiods,1),1);
                for r = 1:size(ripperiods,1)
                    temp = lookup2(ripperiods(r,:),lfpTime);
                    pos = [0; diff(currentDeg(temp(1):temp(2)))];
                    pos(pos<-350) = 0;
                    pos = cumsum(pos);
                    time = sum(diff(vrTime(temp(1):temp(2))));
                    meanspeed = (pos(end)-pos(1)) ./ time;
                    excl6(r) = any(meanspeed > speedTh);
                end
                incl6 = ~excl6;
            end
            
            
            excl = excl1 | excl3 | excl4 | excl5 | excl6;%excl2 removed with tbdratio
            incl = incl1 & incl3 & incl4 & incl5 & incl6;%incl2 removed with tbdratio
            
            
            %find outlier ripples and their index
            if any(excl) && (exclude) %exclude outlier ripples
                rip = ripples(fnum);
                n = fieldnames(rip);
                for i =  1:length(n) %for each field, select only included ripples
                    try
                        eval(['ripples.' n{i} ' = rip.' n{i} '(incl);']); %added semicolon inside ALP 2/14/2020
                    end
                end
            end
            %make list of outlier ripples and index
            excluded = [excluded; excl];
            
            if (exclude) && ~isempty(applyspeed)
                ripples(fnum).excluderipples = ['excluded ', num2str(sum(excl)),...
                    ' ripples with power ratio (', num2str(freqnumerator), '/' ,...
                    num2str(freqdenominator),') < ' num2str(ratiothresh), ' in ', num2str(timearoundrip), ' sec around rip', ...
                    ' and mean speed > ' num2str(speedTh) ' deg/s'];
            else
                ripples(fnum).excluderipples = ['excluded ', num2str(sum(excl)), ...
                    ' ripples with power ratio (', num2str(freqnumerator), '/' , ...
                    num2str(freqdenominator),') < ' num2str(ratiothresh), ' in ', num2str(timearoundrip), ' sec around rip'];
            end
        end
        
        clear excl incl;
    end
clear partindex timearoundrip freqnumerator freqdenominator ratiothresh;
end
