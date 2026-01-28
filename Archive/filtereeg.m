function eeg = filtereeg(datastruct, data, freqband, eegendsamprate)
% filter raw eeg to eeg
% eeg = filtereeg(datastruct, data, index, datadir,  freqband, eegendsamprate)
% filters data for freqband(1)-(2) using an FIR elliptical filter and saves
% as eeg in datadir
% INPUTS
% 	datastruct -- datastruct has info want in results eeg structure
%   data -- data to filter, is Nx1 or 1xN
%   index -- [Animal YYMMDD file#]
%   datadir -- folder to save eeg to
%   freqband -- [1 300]
%
% ASinger 8/24/15

%filter for LFP
% create FIR filter coefficients 1-300Hz
[lfp_b, lfp_a] = ellip(2, 0.1, 40, [freqband]*2/datastruct.samprate);
%for filtering but no decimating: FIR using elliptical impulse
lfp_filtered = filtfilt(lfp_b, lfp_a, data);

%downsample
if datastruct.samprate ~= eegendsamprate
    nsample = datastruct.samprate/eegendsamprate;
    if rem(nsample,1)~=0
        warning('desired downsampled sample rate is not a multiple of original sample rate')
    end
    % down_lfp_filtered = downsample(lfp_filtered,nsample);
    down_lfp_filtered = decimate(lfp_filtered, nsample);
    eeg.downsample = ['downsampled via decimate by ', num2str(nsample)];
else
    down_lfp_filtered = lfp_filtered;
end
%make structure
eeg = datastruct;
eeg.data = down_lfp_filtered;
eeg.filteringinfo = ['bandpass filtered ', num2str(freqband(1)), '-', num2str(freqband(2)), 'Hz      '; 'using an FIR elliptical filter '];
eeg.descript = 'eeg data'; %added due to error for LFP filter
eeg.starttime = 0;
%save structure
%save([datadir, 'eeg', num2str(index(3))], 'eeg')




