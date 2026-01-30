function ripplepostfileprocess2( directoryname, sessindex, files,...
    timearoundrip, freqnumerator, freqdenominator, ratiothresh, outlierindices_allchan, varargin)
% excluded = ripplepostfileprocess2( directoryname, sessindex, timearoundrip, freqnumerator, freqdenominator, ratiothresh, varargin)
% For files with 3 Number index instead of 4.
%inputs
%   ripples - actual ripples
%   lfp data - filtered (eeg)
%   thetas - theta periods
%   tdbratio - ratio theta delta bete
%   need to go back and create all of these
% options -
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

% create the list of files that we should filter
if any(files == -1)
    tmpflist = dir(sprintf('%s*ripples*.mat', directoryname));
    flist = cell(size(tmpflist));
    for i = 1:length(tmpflist)
        flist{i} = sprintf('%s%s', directoryname, tmpflist(i).name);
    end
else
    for f = 1:length(files)
        tmpflist = dir([directoryname, '*ripples',num2str(files(f)),'.mat']);
        if ~isempty(tmpflist)
            flist{f} =  sprintf('%s%s', directoryname, tmpflist(1).name);
        end
    end
end

excluded = [];
% go through each file in flist and filter it
if ~isempty(flist)
    for fnum = 1:length(flist)
        %load files
        load(flist{fnum}); %the ripple file
        filenumber = size(ripples{sessindex(1)}{sessindex(2)},2); %append filenumber
        %         tdbflist = dir(sprintf('%sEEG/tdbratio*%02d.mat', directoryname, filenumber));
        load([directoryname, 'EEG/tdbratio', num2str(filenumber), '.mat']); %load tdbratio file
        index = [sessindex filenumber];
        sprintf(['ripple post processing index: ', num2str(index)]);
        load([directoryname, 'eeg', num2str(index(3))]); %load corresponding eegfile
        load([directoryname, 'thetas', num2str(index(3))]); %load corresponding thetas
        
        if any(ripples{index(1)}{index(2)}{index(3)}.midind)
            %check all samprates =
            if ~isequal(ripples{index(1)}{index(2)}{index(3)}.samprate, thetas{index(1)}{index(2)}{index(3)}.samprate) | ...
                    ~isequal(ripples{index(1)}{index(2)}{index(3)}.samprate, tdbratio{index(1)}{index(2)}{index(3)}.samprate)
                error('sampling rates are not the same')
            end
            
            %get ratio of freq power
            [ratio] = findpowerratioripplevsabove2(index, ripples, eeg, timearoundrip, freqnumerator, freqdenominator);
            if length(ratio) ~= length(ripples{index(1)}{index(2)}{index(3)}.midind)
                error('ratio is not being calculated properly in ripple post processing')
            end
            excl1 = ratio<ratiothresh; %ripples that should b excluded
            incl1 = ratio>=ratiothresh; %ripples that should be included
            
            %find if meet tdbratio criteria as well
            rips = ripples{index(1)}{index(2)}{index(3)}.midind;
            ripperiods = [rips-0.5*ripples{index(1)}{index(2)}{index(3)}.samprate rips+0.5*ripples{index(1)}{index(2)}{index(3)}.samprate]; %take time before and after ripplem so 1 seconds todal
            ripperiods(ripperiods(:,1)<1,1) = 1; %if start before start of recording, set to 1
            ripperiods(ripperiods(:,2)>length(tdbratio{index(1)}{index(2)}{index(3)}.data), 2) = length(tdbratio{index(1)}{index(2)}{index(3)}.data); %if start after end of recording set to end of recording
            
            for r = 1:size(ripperiods)
                meantdb(r,:) = mean(tdbratio{index(1)}{index(2)}{index(3)}.data(ripperiods(r,1):ripperiods(r,2)));%average tdbratio during 1 sec around center of rip
            end
            %             excl2 = meantdb>thetas{index(1)}{index(2)}{index(3)}.baselinethreshold; %exclude if average thetadeltabeta ratio is > thetas.baseline threshold
            %             incl2 = meantdb<=thetas{index(1)}{index(2)}{index(3)}.baselinethreshold;
            %%%% CHANGED TO BASELINE NOT BASELINE THRESH - 6.19.18 SP
            excl2 = meantdb>thetas{index(1)}{index(2)}{index(3)}.baseline; %exclude if average thetadeltabeta ratio is > thetas.baseline 
            incl2 = meantdb<=thetas{index(1)}{index(2)}{index(3)}.baseline;
            
            %excl if outlier falls into ripple period
            %             tempout = getbinindex(eeg{end}{end}{end}.outlierindices, [ripples{index(1)}{index(2)}{index(3)}.startind ripples{index(1)}{index(2)}{index(3)}.endind]);
            tempout = getbinindex(outlierindices_allchan{fnum}, [ripples{index(1)}{index(2)}{index(3)}.startind ripples{index(1)}{index(2)}{index(3)}.endind]);
            excl3 = zeros(size(rips,1),1);
            excl3(unique(tempout(tempout~=0))) = 1;
            incl3 = ~excl3;
            
            %excl if max value is greater than maxvalthresh (default 1500)
            maxvalthresh = 1500;
            tempout2 = zeros(size(rips,1),1);
            for i = 1:length(ripples{index(1)}{index(2)}{index(3)}.startind)
                eegvalues = eeg{index(1)}{index(2)}{index(3)}.data((ripples{index(1)}{index(2)}{index(3)}.startind(i):ripples{index(1)}{index(2)}{index(3)}.endind(i)));
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
                load([directoryname, 'spikes', num2str(index(3))]); %load corresponding eegfile
                spikeidx = spikes{index(1)}{index(2)}{index(3)}.spikeindex;
                downsamp = eeg{index(1)}{index(2)}{index(3)}.downsample;
                rips = [ripples{index(1)}{index(2)}{index(3)}.startind, ripples{index(1)}{index(2)}{index(3)}.endind];
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
                idcs   = strfind(directoryname,filesep);
                rawposdir = directoryname(1:idcs(end-1)-1);
                tempD =  dir(fullfile(rawposdir, '*rawpos*.mat'));
                if isempty(tempD)
                    rawposdir = directoryname(1:idcs(end-2)-1); %NJ added 05/17/22 bc rawpos files may be outside brainReg folder
                end
                rawpos = load(fullfile(rawposdir, ['rawpos' num2str(index(3)) '.mat']));
                rawpos = rawpos.rawpos{index(1)}{index(2)}{index(3)};
                downsamp = eeg{index(1)}{index(2)}{index(3)}.downsample;
                rips = ripples{index(1)}{index(2)}{index(3)}.midind;
                rips = rips.*downsamp; %turn it into spike index at 30000 samprate
                ripperiods = [rips-t_sec*ripples{index(1)}{index(2)}{index(3)}.samprate*downsamp, rips+t_sec*ripples{index(1)}{index(2)}{index(3)}.samprate*downsamp]; %take time before and after ripplem so 2 seconds total - had to make it 30k samp rate bc using spike index
                ripperiods(ripperiods(:,1)<1,1) = 1; %if start before start of recording, set to 1
                ripperiods(ripperiods(:,2)>length(tdbratio{index(1)}{index(2)}{index(3)}.data)*downsamp, 2) = length(tdbratio{index(1)}{index(2)}{index(3)}.data)*downsamp; %if start after end of recording set to end of recording
                
                excl6 = false(size(ripperiods,1),1);
                for r = 1:size(ripperiods,1)
                    temp = lookup2(ripperiods(r,:),rawpos.ephysInd);
                    pos = [0; diff(rawpos.currentTheta(temp(1):temp(2)))];
                    pos(pos<-350) = 0;
                    pos = cumsum(pos);
                    time = sum(diff(rawpos.vrTime(temp(1):temp(2))));
                    meanspeed = (pos(end)-pos(1)) ./ time;
                    excl6(r) = any(meanspeed > speedTh);
                end
                incl6 = ~excl6;
            end
            
            
            excl = excl1 | excl2 | excl3 | excl4 | excl5 | excl6;
            incl = incl1 & incl2 & incl3 & incl4 & incl5 & incl6;
            
            
            %find outlier ripples and their index
            if any(excl) && (exclude) %exclude outlier ripples
                rip = ripples{index(1)}{index(2)}{index(3)};
                n = fieldnames(rip);
                for i =  1:length(n) %for each field, select only included ripples
                    try
                        eval(['ripples{index(1)}{index(2)}{index(3)}.' n{i} ' = rip.' n{i} '(incl);']); %added semicolon inside ALP 2/14/2020
                    end
                end
            end
            %make list of outlier ripples and index
            excluded = [excluded; repmat(index, length(excl), 1) excl];
            
            if (exclude) && ~isempty(applyspeed)
                ripples{index(1)}{index(2)}{index(3)}.excluderipples = ['excluded ', num2str(sum(excl)),...
                    ' ripples with power ratio (', num2str(freqnumerator), '/' ,...
                    num2str(freqdenominator),') < ' num2str(ratiothresh), ' in ', num2str(timearoundrip), ' sec around rip', ...
                    ' and mean speed > ' num2str(speedTh) ' deg/s'];
            else
                ripples{index(1)}{index(2)}{index(3)}.excluderipples = ['excluded ', num2str(sum(excl)), ...
                    ' ripples with power ratio (', num2str(freqnumerator), '/' , ...
                    num2str(freqdenominator),') < ' num2str(ratiothresh), ' in ', num2str(timearoundrip), ' sec around rip'];
            end
            %save new ripple structure
            save(flist{fnum}, 'ripples');
        end
        
        clear ripples thetas tdbratio meantdb ratio excl incl eeg;
    end
end
clear directoryname partindex timearoundrip freqnumerator freqdenominator ratiothresh;
end
