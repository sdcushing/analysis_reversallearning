function [dirs, params] = getDirectoriesAndParams_JLK_DC()
%% JLK default directories for analyses, last checked 10/15/24
%updated from NJ
%ephys
dirs.rawdata = '\\ad.gatech.edu\bme\labs\singer\RawData\DAlesion\';
dirs.processeddata = '\\ad.gatech.edu\bme\labs\singer\ProcessedData\VR_Novelty_DAlesion\';
dirs.filters = '\\ad.gatech.edu\bme\labs\singer\Data Extraction and PreProcessing\Intan\Filtering\';
dirs.code = '\\ad.gatech.edu\bme\labs\singer\Danielle\code\AnalysisCode\Neuropixels_analyses\';
%behavior
dirs.virmenrawdata = '\\ad.gatech.edu\bme\labs\singer\Danielle\Behavior\Annular FAM\';
dirs.virmensessiondata = '\\ad.gatech.edu\bme\labs\singer\Danielle\Behavior\sessionData\';
dirs.rewardzones = '\\ad.gatech.edu\bme\labs\singer\Danielle\Behavior\TrackFigures\RewardZones.mat';
%save results
dirs.savefigures = '\\ad.gatech.edu\bme\labs\singer\Danielle\Figures\';
dirs.savefiguresbeh = '\\ad.gatech.edu\bme\labs\singer\Danielle\Figures\behavior\';
dirs.saveoutputstructs = '\\ad.gatech.edu\bme\labs\singer\Danielle\OutputStructs\';
%spreadsheets with training or recording info
%dirs.spreadsheet_rec = '\\ad.gatech.edu\bme\labs\singer\Josh\Spreadsheets\UpdateBehaviorSpreadsheet.csv'; %JLK commented out 10/2/24 because currently only one spreadsheet
dirs.spreadsheet = '\\ad.gatech.edu\bme\labs\singer\Danielle\Behavior\Annular FAM\NoveltyBehaviorSpreadsheet2.xlsx'; 
dirs.clusfolder = 'sorted\';
dirs.cluster_local = 'C:\Users\scushing6\Desktop\TempKilosort\';
dirs.kilosortPyEnv = 'C:\Users\scushing6\AppData\Local\anaconda3\envs\kilosort\python.exe';
dirs.kilosortPyScript = 'Y:\singer\Danielle\Code\AnalysisCode\Neuropixels_analyses\runKilosort4.py';

%% parameters for analyses
%general
params.iden = 'DC'; %JK for default Josh's mice
params.animals = [21];
params.datesincl = [];
params.datesexcl = [];
params.recday = [];
params.brainReg = {'CA1'};
params.probeChannels = {1:64};                                              %64-channel NeuroNexus probe
params.binsize_ms = 1;                                                      %in ms for monoconnex calculations
params.binsize_deg = 5;                                                     %in degrees for behavioral analyses
params.binsize_s = 60;                                                      %in sec for firing rate stability across session
params.samprate = 20000;  
params.lfp_samprate_down = 2000;%in Hz for SpikeGadgets acquisition system
params.timeoutZone = 17;
params.rocID = {'speed', 'lickrate', 'deltalickrate'};
params.rocMultiplier = [-1 1 1];
params.numTrPerBlock = 5;
params.savechnum = {0:63}; %how to save the channels listed below in channelInds. Will have portRegion as external folder

%settings from og preprocessing pipeline to get rhd2mat to work
params.ttlCh = [1 2 3];
params.triggerfilename = 'sync';
params.ttlfilename = {'reward'};
params.autoTriggerCh = 4;
params.recordingTriggerCh = 0;
params.exportbehaviorfunc = 'rawpos_default_DC';
params.adc_channels = [0 1];
params.adcInfo = 'scaled ball tracking';
params.velocityfilename = 'velocity';
params.savechnum = {0:63};
params.rm60HzNoise = 0;

%zones
params.cueSize = 10;                                                        %in degrees, size of wall cue (zone)
params.gapBefore = 3 * params.cueSize;
params.gapAfter = 3 * params.cueSize;

%set colors for original and update sessions
up_colors = cbrewer('seq','Greens',12);
params.colors_nov = up_colors(6:2:end, :); %shades of green for update
og_colors = cbrewer('seq','Greys',12);
params.colors_fam = og_colors(7:2:end, :); %shades of grey for original
% goal_colors = cbrewer('seq','Blues',8);
% params.colors_goalstim = goal_colors(3:2:end, :); %shades of blue for goal stim 
% sham_colors = cbrewer('seq','Oranges',8);
% params.colors_shamstim = sham_colors(3:2:end, :); %shades of orange for sham stim
% JLK commented out lines below 10/2/24
% params.colors_shamstim = hex2rgb({'#ffba61','#ffa42e','#fa8d00'});
% params.colors_goalstim = hex2rgb({'#61a5ff','#2e88ff','#006cfa'});
% params.colors_narrowInt = hex2rgb({'#d2dbec','#7894c5','#36489b'});
% params.colors_wideInt = hex2rgb({'#e3f4f7','#aadde8','#72c7d9'});
% params.colors_pyr = hex2rgb({'#f5d6d6','#e08585','#cb3433'});

%ripples
params.outliernstd = 15; %original: 15, strict: 7 - number of std to exclude outliers. [1x1] or [1x2]. see interpoveroutliers for details
params.filteegfreq = [1 300];
params.eegsamprate = 2000;
params.ripple.nstdEnv = [3 3];%number of standard deviations above mean of envelope for 150-250Hz to detect ripple
params.ripple.nstdNoise = 10;%detect noise deflections in raw data
params.ripple.minRipDur = 0.02;%ms
params.ripple.timeAroundRip = 0.05;%time (ms) around midpoint of rip to compute ratio
params.ripple.freqNumerator = [100 250];%high frequency numerator
params.ripple.freqDenominator = [250 400];%high frequency denominator
params.ripple.ratioThresh = 2;%lowered from 5
params.ripple.applyCriteria = 1;
params.ripple.applyHighFreqRatio = 0;
params.ripple.applyTDBRatio = 0;
params.ripple.applySpeed = 1;
params.ripple.applyMUA = 1;

%cell type identification
params.spikewidthTh = 0.5; %in ms
params.autocorrTh = 4.5; %in ms
params.intRateTh = 5; %in Hz
params.pyrRateTh = 10; %in Hz
params.snrTh = 1; %signal-to-noise threshold 
params.ttlsamprate = 10000;

%place coding
params.environments = {'og','up', 'nov', 'nov2'};
params.speedTh = 2;                                                         %in deg/s for occupancy normalized rate maps
params.spatialinfoTh = 95;                                                  %in n'th percentile for sptial info
params.plotCellT = 'Interneuron';                                           %default cell-type to plot 

%parameters for monoconnex
params.nonnorm = 0; %raw xcorr
params.corrcoeffnorm = 0; %xcorr normalized by norm(spiketraincella)*norm(spiketraincellb)
params.geommeanfrnorm = 1; %xcorr normalized by sqrt(meanfrcella*meanfrcellb)
params.rawpeakvals = 0;
params.efficacycontribtionpeakvals = 0;
params.meanbaselinepeakvals = 1;
params.integralbaselinepeakvals = 0;
params.rangeforbaseline = 2; %how many bins to use to measure baseline, want the range*2 to be the same number you integrate over for the peak part
params.wheretomeasurebaseline = 1; %how far +/- the centerbin:centerbin+rangeforbaseline*2 area (where you meausre the peak) that you want to measure baseline from
%changed this from 1
params.rectify = 1;
params.lagtime = 50;
params.longtimescales = 0;
params.refactoredcalc = 0;

%to overwrite existing files or not
params.rewrite.behavior = 0;
params.rewrite.ROC = 1;
params.rewrite.clusters = 0;
params.rewrite.cell_metrics = 0;
params.rewrite.neuralStructs = 1;
params.rewrite.cellYield = 0;
params.rewrite.pyrLayer = 0;
params.rewrite.ripples = 1;
params.rewrite.decodingData = 0;
params.rewrite.decodingResults = 0;
