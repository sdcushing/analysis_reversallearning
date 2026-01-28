function [phase, env, filtdata] = filtereeg2_DC(eeg, kernel)
%based on filtereeg2_intan, but compatible with getRipplesDC
%inputs
%   eeg - filtered lfp data, size channelsxsamples
%   kernal - kernal for the filter
for chCtr = 1:size(eeg, 1)
    %ripple filter (150-250 Hz)
    filtdata(chCtr,:) = filtfilt(kernel, 1 , eeg(chCtr,:));

    %Hilbert transform ripple data
    hdata = hilbert(filtdata(chCtr,:));
    phase(chCtr,:) = angle(hdata);
    env(chCtr,:) = abs(hdata);
    clear hdata
end

end