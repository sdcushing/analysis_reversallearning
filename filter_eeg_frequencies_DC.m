function [rawDataBySessionNeural] = filter_eeg_frequencies_DC(rawDataBySessionNeural, dirs, params, saveNeuralPath)
%adapted from filtereeg2_Intan.m, extractripples3.m,
% ripplepostfileprocess2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% initialize params and load filters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = [];
beta = [];
delta = [];
tdbratio = [];
ripple = [];
chCtr = 0;
layerChans = sessionPyrLayerInfo.pyrLayerCA1;%CA1 only
%layerChans = [sessionPyrLayerInfo.pyrLayerCA3 sessionPyrLayerInfo.pyrLayerCA1];%CA3+CA1
load([dirs.code 'ripplefilter.mat'])
load([dirs.code 'thetafilter.mat'])
load([dirs.code 'betafilter.mat'])
load([dirs.code 'deltafilter.mat'])

%%%%% add filter description and parameters %%%%%
ripple.descript = ripplefilter.descript;
ripple.kernel = ripplefilter.kernel;
ripple.samprate = ripplefilter.samprate;
smoothing_width_rip = 0.004; % 4 ms
smoothing_kernel_rip = gaussian(smoothing_width_rip*ripple.samprate, ceil(8*smoothing_width_rip*ripple.samprate));
smoothing_width_tdb = 1; % 1 s
smoothing_kernel_tdb = gaussian(smoothing_width_tdb*ripple.samprate, ceil(8*smoothing_width_tdb*ripple.samprate));
minRipDur = round(params.ripple.minRipDur * ripple.samprate);
tmpDat = rawDataBySessionNeural.lfpData;

%josh has some noise removal here. not going to worry about that right now,
%but should consider at some point. I have outlier removal when i filter +
%downsample lfp in getneuralstructs

%make ripple filtered eeg traces - based on ripplefileprocess2 and
%filtereeg2_intan
if ~(isfield(rawDataBySessionNeural, 'ripple'))
[ripple.phase, ripple.env, ripple.data, ripple.samprate] = filtereeg2_DC(tmpDat, rawDataBySessionNeural.lfp_meta.downsamplerate, ripplefilter);
%getting other frequencies for ripple outlier removal
[theta.phase, theta.env, theta.data, theta.samprate] = filtereeg2_DC(tmpDat, rawDataBySessionNeural.lfp_meta.downsamplerate, thetafilter);
[beta.phase, beta.env, beta.data, beta.samprate] = filtereeg2_DC(tmpDat, rawDataBySessionNeural.lfp_meta.downsamplerate, betafilter);
[delta.phase, delta.env, delta.data, delta.samprate] = filtereeg2_DC(tmpDat, rawDataBySessionNeural.lfp_meta.downsamplerate, deltafilter);
tdbratio.data = theta.data./(beta.data + delta.data);
    smoothing_width = 1; % % define the standard deviation for the Gaussian smoother
    kernel = gaussian(smoothing_width*rawDataBySessionNeural.lfp_meta.downsamplerate, ...
        ceil(8*smoothing_width*rawDataBySessionNeural.lfp_meta.downsamplerate));
for i = 1:size(tdbratio, 1)
    tdbratio.data(i, :) = smoothvect(tdbratio.data(i,:), kernel);
end
rawDataBySessionNeural.ripple = ripple;
rawDataBySessionNeural.theta = theta;
rawDataBySessionNeural.beta = beta;
rawDataBySessionNeural.delta = delta;
rawDataBySessionNeural.tdbratio = tdbratio;
save([saveNeuralPath '\rawDataBySessionNeural.mat'], 'rawDataBySessionNeural', '-v7.3');
else
    ripple = rawDataBySessionNeural.ripple;
    theta = rawDataBySessionNeural.theta;
    beta = rawDataBySessionNeural.beta;
    delta = rawDataBySessionNeural.delta;
    tdbratio = rawDataBySessionNeural.tdbratio;
end

end