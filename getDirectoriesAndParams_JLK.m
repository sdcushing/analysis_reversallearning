function [dirs, params] = getDirectoriesAndParams_JLK()
%% JLK default directories for analyses, last checked 10/15/24
%updated from NJ
%ephys
dirs.rawdata = '\\ad.gatech.edu\bme\labs\singer\RawData\VR_AnnularReversalTask_JK\';
dirs.processeddata = '\\ad.gatech.edu\bme\labs\singer\ProcessedData\VR_AnnularReversalTask_JK\';
dirs.code = '\\ad.gatech.edu\bme\labs\singer\Josh\Code\AnalysisCode\Neuropixels_analyses\';
dirs.kilosortPyEnv = 'C:\Users\jkrasney3\AppData\Local\miniconda3\envs\kilosort\python.exe';
dirs.kilosortPyScript = 'Y:\singer\Josh\Code\AnalysisCode\Neuropixels_analyses\runKilosort4.py'; 
%behavior
dirs.virmenrawdata = '\\ad.gatech.edu\bme\labs\singer\Josh\Behavior\AnnularFAM\';
dirs.rewardzones = '\\ad.gatech.edu\bme\labs\singer\Josh\Behavior\TrackFigures\RewardZones.mat';
%save results
dirs.savefigures = '\\ad.gatech.edu\bme\labs\singer\Josh\Figures\';
dirs.saveoutputstructs = '\\ad.gatech.edu\bme\labs\singer\Josh\OutputStructs\';
%spreadsheets with training or recording info
dirs.spreadsheet = '\\ad.gatech.edu\bme\labs\singer\Josh\Spreadsheets\UpdateBehaviorSpreadsheet.csv'; 
dirs.clusfolder = 'sorted\';
dirs.cluster_local = 'C:\Users\njeong9\Documents\TEST\';

%% parameters for analyses
%general
params.iden = 'JK'; %JK for default Josh's mice
%params.animals = [5 6 7 9 10 12 13 14 15 18 19 20 21 22 23 24 25];%animals with usable VR data
%params.animals = [14 15 18 19 20 21 22 23 24 25];%tested animals
params.animals = [14];
params.controlGroup = [7 8 12 20 21 22 23 24 25];%APOE3
params.experimentalGroup = [5 6 9 10 12 13 14 15 18 19];%APOE4
params.testedMice = [14 15 18 19 20 21 22 23 24 25];
params.datesincl = [];
params.datesexcl = [];
params.recday = [];
params.brainReg = {'CA1'};
params.probeChannels = {1:64};                                              %64-channel NeuroNexus probe
params.binsize_ms = 1;                                                      %in ms for monoconnex calculations
params.binsize_deg = 2;                                                     %in degrees for behavioral analyses
params.binsize_stime = 6;                                                   %in sec for neural analyses
params.binsize_mstime = 50;                                                  %in sec for neural analyses
params.binsize_s = 60;                                                      %in sec for firing rate stability across session
params.samprate = 30000;                                                    %in Hz for SpikeGadgets acquisition system
params.timeoutZone = 17;
params.rocID = {'speed', 'lickrate', 'deltalickrate'};
params.rocMultiplier = [-1 1 1];
params.numTrPerBlock = 5;

%zones
params.OGrZones = [30,140,240];%original reward zones
params.UPrZones = [90,200,280];%current updated reward zones
params.cueSize = 10;                                                        %in degrees, size of wall cue (zone)
params.gapBefore = 2 * params.cueSize;
params.gapAfter = 2 * params.cueSize;

%neural
params.ap_samprate = 30000;
params.lfp_samprate = 2500;
params.lfp_samprate_down = 2000;
params.rm60HzNoise = 0;
params.nstdNoise = 6;%JLK changed from 4 12/4/25 after comparing noise across channels
params.numChans = 383;%0-based
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

%place coding
params.environments = {'og','up', 'nov', 'nov2'};
params.speedTh = .3;                                                         %in deg/s
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
params.rewrite.decodingData = 1;
params.rewrite.decodingResults = 1;
