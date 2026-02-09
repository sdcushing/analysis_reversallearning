function plotLFPperiods(savefigsdir, plottingdatadir, sessindex, files, params, ...
    plotEEG, plotthetas, plotnonthetas, plotripples, plotgammas, brainReg, ...
    bestRippleChan, interactive)
%plotLFPperiods
% ALP 4/01/2021
% this function plots LFP periods of interest for manual inspection of the
% data. 

plotthetasch = bestRippleChan; 
plotnonthetasch = bestRippleChan; 
plotripplesch = bestRippleChan; 
plotgammasch = bestRippleChan; 

if plotEEG
    for f = 1:length(files)
        index = [sessindex files(f)]; 
        ploteegartifacts_041819(savefigsdir, plottingdatadir, index, brainReg, interactive)
    end
    close all
end

if plotthetas
    %removed plotting comparing theta thresholds, to get that code see script pre-20200108
    plotexamples_thetaperiods_gamma_extracell(savefigsdir, plottingdatadir, sessindex, ...
        files, 'thetas','theta', 'lowgamma', 1.5, plotthetasch, brainReg, interactive)
    close all
end

if plotnonthetas
    plotexamples_thetaperiods_extracell(savefigsdir, plottingdatadir, sessindex, files, ...
        'nonthetas', 3, plotnonthetasch, brainReg, interactive)
    close all
end

if plotripples
    timearoundrip = 1; %for plotting
    [mnP, maxP, ~] = plotexamples_ripplesenergygram2(savefigsdir, plottingdatadir, ...
        sessindex, files, timearoundrip, params.extractripples_freqnumerator, ...
        params.extractripples_freqdenominator,[80 450], 1, brainReg, interactive);
    
    ratiolabel = [num2str(params.extractripples_freqnumerator(1)), '-', ...
        num2str(params.extractripples_freqnumerator(2)), '/', ...
        num2str(params.extractripples_freqdenominator(1)), '-', ...
        num2str(params.extractripples_freqdenominator(2)) ] ;
    xlabel(['mean Power ratio: ', ratiolabel, '100ms centered at center of ripple'])
    close all
end

if plotgammas
    gammasubsamp = 300; timearoundgamma = 1; %for plotting
    plotexamples_gammas(savefigsdir, plottingdatadir, sessindex, files, timearoundgamma, ...
        gammasubsamp, params.extractgamma_altnstd, brainReg, interactive)
    close all
    
end


end

