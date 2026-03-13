function [rawDataBySessionNeural] = getRipples_DC3(dirs, params, saveNeuralPath, plotRipples, subj, sessNum)
%adapted from filtereeg2_Intan.m, extractripples3.m,
% ripplepostfileprocess2, findpowerratioripplevsabove2
%last checked JLK 1/8/26
%DC changing to fully incorporate ripplefileprocess, extractripples3, 
%ripplepostfileprocess2, getBestRippleChan_simple, plotLFPperiods


%% extract ripples across session %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load session data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%need to include ripplesettings.stdev, ripplesettings.base to pass to
%%extractripples3
load([saveNeuralPath '\rawDataBySessionNeural.mat'])
load([saveNeuralPath '\sessionPyrLayerInfo.mat'])
load([saveNeuralPath '\ripplesettings.mat'])

%removing old ripples
if isfield(rawDataBySessionNeural,'ripplesBad')
    rawDataBySessionNeural = rmfield(rawDataBySessionNeural,'ripplesBad');
end
if isfield(rawDataBySessionNeural,'ripplesGood')
    rawDataBySessionNeural = rmfield(rawDataBySessionNeural,'ripplesGood');
end

    ripple = rawDataBySessionNeural.ripple;
    theta = rawDataBySessionNeural.theta;
    beta = rawDataBySessionNeural.beta;
    delta = rawDataBySessionNeural.delta;
    tdbratio = rawDataBySessionNeural.tdbratio;

%actually detect ripples - updated extractripples3 rather then
%retype it all
ripples = extractripples3_DC(ripple, params.ripple.minRipDur, params.ripple.nstdEnv, ...
    'stdev', ripplesettings.stdev, 'baseline', ripplesettings.base);%options for more inputs, but not fully
%functional yet. on DC todo list

%bad ripples are ripples that do not pass criteria. based on
%outlierindices_allchan and ripplepostfilesprecess2
%find outliers
outlierindices_allchan = rip_outlierexcludeallchan_DC(rawDataBySessionNeural.lfp_meta);%fill out

%based on ripplepostfileprocess2, updated to be compatible with this
%pipeline. still need to add pos and mua functionality. need to make sure
%output is consistent with JK output (badRipples and goodRipples)
%for pos (speed) will need currentDeg and lfpTime (already in correct
%sample rate) then look up temp = lookup2(ripperiods(r,:),rawpos.ephysInd);
%and find pos at those indices, get speed from that, if > threshold exclude
ripples_bad = [];
ripplepostfileprocess2_DC(ripples, tdbratio, tmpDat, params.ripple.timeAroundRip, params.ripple.freqNumerator, ...
    params.ripple.freqDenominator, params.ripple.ratioThresh, outlierindices_allchan, ...
    rawDataBySessionNeural.currentDeg, rawDataBySessionNeural.vrTime, rawDataBySessionNeural.lfpTime, ripples_bad, 'exclude', 1, 'applyspeed', params.rippostprocess_applySpeed);
%based on getbestripplechannelsimple. planning on plotting ~raw trace
%(downsampled with outlier filter) as well as ripple itself
if plotRipples
    ripple.bestchan = getBestRippleChan_DC(ripple, ripples);
    plottingchan = ripple.bestchan.channel;
    ripple.bestchan.channum1indx = rawDataBySessionNeural.lfp_meta.channelInd(plottingchan);
    %plot everything of interest on the best ripple channel
    savefigsdir = fullfile(dirs.savefigures, 'ProcessingFigures', params.iden, params.brainReg, 'ripples', sessNum);
    %plottingdatadir = [probeprocesseddatadir, num2str(plottingChan), '\'];  
%%%%%%%%%%%
    timearoundrip = 1; plotexamples = 1;%for plotting
    interactive = 0;
    [mnP, maxP] = plotexamples_ripplesenergygram2_DC(savefigsdir, ripple, ripples, tmpDat, timearoundrip, params.ripple.freqNumerator,...
    params.ripple.freqDenominator, [80 450], interactive, subj, plotexamples, rawDataBySessionNeural.speedSmooth, rawDataBySessionNeural.lfpTime, ...
    params.rippostprocess_applySpeed);

    ratiolabel = [num2str(params.ripple.freqNumerator(1)), '-', ...
        num2str(params.ripple.freqNumerator(2)), '/', ...
        num2str(params.ripple.freqDenominator(1)), '-', ...
        num2str(params.ripple.freqDenominator(2)) ] ;
    xlabel(['mean Power ratio: ', ratiolabel, '100ms centered at center of ripple'])
    close all
end
%add all the things to the rawneural structure and save
rawDataBySessionNeural.ripples = ripples;
rawDataBySessionNeural.ripple = ripple;
rawDataBySessionNeural.ripplesBad = ripples_bad;
save([saveNeuralPath '\rawDataBySessionNeural.mat'], 'rawDataBySessionNeural', '-v7.3');
end

