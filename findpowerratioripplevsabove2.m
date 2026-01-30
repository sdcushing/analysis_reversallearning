function [ratio] = findpowerratioripplevsabove2(ripples, eeg, timearoundrip, freqnumerator, freqdenominator)
% For files with 3 Number index instead of 4.

%make times vector
%ALP 3/31/2021 edited to be max to account for column or row vectors
samprate = ripples.samprate;

a = max(size(eeg.data))/samprate;
times = [0:1/samprate:a]';
times = times(1:end-1);

%get ripple indices + time around
rips = ripples.midind;
ripperiods = [rips-timearoundrip*samprate rips+timearoundrip*samprate]; %take time before and after ripplem so ~2 seconds todal

mnP = [];
while ~isempty(ripperiods) & ripperiods(1,1)<0 %if first ripple too close to start, exclude it
    %ripperiods = ripperiods(2:end,:);
    ripperiods(1,:) = [nan nan];
end
while ~isempty(ripperiods)& ripperiods(end,2) > size(times,1) %if last ripple too close to end, exclude it
    %ripperiods = ripperiods(1:end-1,:);
    ripperiods(end,:) = [nan nan];
end

%for each SWR
for r = 1:size(ripperiods,1)
    ripind = [ripperiods(r,1): ripperiods(r,2)]; %indices to select for rip, includes time around ripplss and combine overlapping periods
    
    if max(ripind)<size(times,1) && min(ripind) > 0 %added  min(ripind) > 0 bc rare case where ripple period starts from 0, NJ 6/16/22
        %compute psd
        [Pxx,F] = pwelch(detrend(eeg.data(ripind)),[],[],[],samprate);
        
        %compute max & mean power in 150-250 and 250-350
        subfreqs{1} = find(F>=freqnumerator(1) & F<=freqnumerator(2)); %ripple band
        subfreqs{2} = find(F>freqdenominator(1) & F<=freqdenominator(2)); %above ripple, prob movement
        for f = 1:2 %for 150-250 or 250-400
            tempmnP(1,f) = mean(Pxx(subfreqs{f}));
            % tempmaxP(1,f) = max(Pxx(subfreqs{f}));
        end
        mnP = [mnP; tempmnP];
        %maxP = [maxP; tempmaxP];
    else
        mnP = [mnP; nan nan];
    end
end
if ~isempty(ripperiods)
    ratio = mnP(:,1)./mnP(:,2);
else
    ratio = [];
end