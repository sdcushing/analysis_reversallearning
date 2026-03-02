function [phase, env, filtdata, samprate] = filtereeg2_DC(eeg, samplerate, filter)
%based on filtereeg2_intan, but compatible with getRipplesDC
%inputs
%   eeg - filtered lfp data, size channelsxsamples
%   kernal - kernal for the filter
if (~isfield(filter, 'samprate'))
    filter.samprate = filter.Fs;
end
if samplerate ~= filter.samprate
    error('filterstruct.samprate does not match eegstruct.samprate');
end
for chCtr = 1:size(eeg, 1)
    %ripple filter (150-250 Hz)
    samprate = samplerate;
    if (isfield(filter, 'kernel'))
        filtdata(chCtr,:) = filtfilt(filter.kernel, 1 , eeg(chCtr,:));
    elseif (isfield(filter, 'tf'))
        filtdata(chCtr,:) = filtfilt(filter.tf.num, filter.tf.den, ...
        eeg(chCtr,:));
    end

    %Hilbert transform ripple data
    hdata = hilbert(filtdata(chCtr,:));
    phase(chCtr,:) = angle(hdata);
    env(chCtr,:) = abs(hdata);
    clear hdata
end

end