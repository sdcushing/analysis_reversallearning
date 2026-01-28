function outlierindices_allchan = rip_outlierexcludeallchan_DC(lfp_meta)
% For files with 3 Number index instead of 4.
% OUTPUTS
%   outlierindices -     cell array of outlier indices for each file
%
% run before ripplepostfileprocess2 to exclude ripples that occur during
% outliers on any channel
%inputs
%   

    %loop through all of the channels for that file
for c = 1:size(lfp_meta.outlierindices, 1) %loop through all channels

        %find  unique outliers for each file
        if any(lfp_meta.outlierindices{i})
            outlierstemp = [outlierstemp, lfp_meta.outlierindices{i}];
        end
    end
    outliers_allchan = unique(outlierstemp);
    outlierindices_allchan = outliers_allchan;
end
