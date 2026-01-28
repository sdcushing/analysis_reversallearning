function outlierindices_allchan = ripplepostfileprocess_outlierexcludeallchan(processeddatadir, sessindex,  files, lfpChans)
% For files with 3 Number index instead of 4.
% OUTPUTS
%   outlierindices -     cell array of outlier indices for each file
%
% run before ripplepostfileprocess2 to exclude ripples that occur during
% outliers on any channel

% loop through all of the files from that recording
for f = 1:length(files)
    ['file num' num2str(f)];
    outlierstemp = [];
    %loop through all of the channels for that file
    %-- updated loop ALP 4/8/21 to account for cases where the channel # might
    %not start at 0 for some reason (for example, some channels disabled,
    %etc) -- 
    for c = 1:length(lfpChans) 
        chanIdx = lfpChans(c); 
        anprocesseddatadir = [processeddatadir, num2str(chanIdx), '/'];
        tmpflist = dir([anprocesseddatadir, '*eeg',num2str(files(f)),'.mat']);
        if ~isempty(tmpflist)
            flist{f} =  sprintf('%s%s', anprocesseddatadir, tmpflist(1).name);
        end

        %load files
        load(flist{f}) %the eeg file
        filenumber = size(eeg{sessindex(1)}{sessindex(2)},2); %append filenumber
        index = [sessindex filenumber chanIdx];

        %find  unique outliers for each file
        if any(eeg{index(1)}{index(2)}{index(3)}.outlierindices)
            outlierstemp = [outlierstemp, eeg{index(1)}{index(2)}{index(3)}.outlierindices];
        end
    end
    outliers_allchan = unique(outlierstemp);
    outlierindices_allchan{f} = outliers_allchan;
end
end
