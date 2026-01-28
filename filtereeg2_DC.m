function filteredeeg = filtereeg2_DC(eeg, kernel, smoothing_kernel_rip)
%based on filtereeg2_intan, but compatible with getRipplesDC
%inputs
%   eeg - filtered lfp data, size channelsxsamples
%   kernal - kernal for the filter
filteredeeg = [];
for chCtr = 1:size(eeg, 1)
    %ripple filter (150-250 Hz)
    filteredeeg.filtdata(chCtr,:) = filtfilt(kernel, 1 , eeg(chCtr,:));

    %Hilbert transform ripple data
    hdata = hilbert(filteredeeg.filtdata(chCtr,:));
    filteredeeg.phase(chCtr,:) = angle(hdata);
    filteredeeg.env(chCtr,:) = abs(hdata);
    filteredeeg.env(chCtr,:) = smoothvect(filteredeeg.env(chCtr,:), smoothing_kernel_rip);%smooth env
    clear hdata
end
end