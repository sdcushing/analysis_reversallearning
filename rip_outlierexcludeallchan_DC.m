function outlierindices_allchan = rip_outlierexcludeallchan_DC(lfp_meta)
% For files with 3 Number index instead of 4.
% OUTPUTS
%   outlierindices -     cell array of outlier indices for each file
%
% run before ripplepostfileprocess2 to exclude ripples that occur during
% outliers on any channel
%inputs
%   
outlierstemp = [];
    %loop through all of the channels for that file
for c = 1:size(lfp_meta.outlierindices, 2) %loop through all channels

        %find  unique outliers for each file
        if any(lfp_meta.outlierindices{c})
            outlierstemp = [outlierstemp, lfp_meta.outlierindices{c}];
        end
        if isnan(outlierstemp)
            outlierstemp = [];
        end

    outlierindices_allchan{c} = outlierstemp;
end

end
