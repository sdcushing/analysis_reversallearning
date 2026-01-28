function [eeg, outlierthresh, numoutlierperiods, outlierindices, outlierperiods] = makefiltereeg_Intan_DC(raweeg, freq, eegsamprate, outliernstd)
%
% INPUTS
% freq -- frequencies for bandpass filtering, usually [1 300] for typical
%   eeg
% eegsamprate -- final sampling rate of eeg, usually 2000
% Updated ALP 4/18/19 to include outlier periods in saving
%updated DC 1/27/26 for compatibility with new pipeline

samprate = raweeg.samprate;

%replace outliers in LFP 
[neweeg, outlierindices, outlierthresh, outlierperiods, numoutlierperiods, ~] = interpoveroutliers_041819(raweeg.data, outliernstd, 2000); %changed from 2 (for mV) to 2000 (for uV)

%filter eeg
eeg = raweeg;
eeg.outlierindices = outlierindices;
eeg.outlierthresh = outlierthresh;
eeg.outlierperiods = outlierperiods;
eeg.numoutlierperiods = numoutlierperiods;


eeg = filtereeg(eeg, neweeg, [freq], eegsamprate);


