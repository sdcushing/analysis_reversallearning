function [rawDataBySessionNeural] = getRipples_DC2(dirs, params, saveNeuralPath, plotRawTraceForRipples, plotRipples)
%adapted from filtereeg2_Intan.m, extractripples3.m,
% ripplepostfileprocess2, findpowerratioripplevsabove2
%last checked JLK 1/8/26
%DC changing to fully incorporate ripplefileprocess, extractripples3, 
%ripplepostfileprocess2, getBestRippleChan_simple, plotLFPperiods

%JK has functionality to do theta etc filtering, should be easy to include
%here too

%% extract ripples across session %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load session data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([saveNeuralPath '\rawDataBySessionNeural.mat'])
load([saveNeuralPath '\sessionPyrLayerInfo.mat'])

%removing old ripples
if isfield(rawDataBySessionNeural,'ripplesBad')
    rawDataBySessionNeural = rmfield(rawDataBySessionNeural,'ripplesBad');
end
if isfield(rawDataBySessionNeural,'ripplesGood')
    rawDataBySessionNeural = rmfield(rawDataBySessionNeural,'ripplesGood');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% initialize params and load filters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = [];
beta = [];
delta = [];
tdbratio = [];
ripple = [];
ripples = [];
ripplesBad = [];
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
%filtereeg2_intan. add functionality for other frequency bands
[ripple.phase, ripple.env, ripple.filtdata] = filtereeg2_DC(tmpDat, ripplefilter.kernel, smoothing_kernel_rip);

%actually detect ripples - updated extractripples3 rather then
%retype it all
ripples = extractripples3_DC(ripple, params.ripple.minRipDur, params.ripple.nstdEnv)%options for more inputs, but not fully
%functional yet. on DC todo list

%bad ripples are ripples that do not pass criteria. based on
%outlierindices_allchan and ripplepostfilesprecess2
%find outliers
outlierindices_allchan = rip_outlierexcludeallchan_DC(rawDataBySessionNeural.lfpmeta);%fill out


end

