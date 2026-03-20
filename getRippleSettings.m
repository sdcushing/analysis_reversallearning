function [ripplesettings] = getRippleSettings(saveNeuralPath, files, sessDate)
%Function to combine all ripple envelopes for all sessions in a given day
%and create the threshold, baseline, and stdev
%   improved performance for detection of ripples in animals who rarely
%   stop
full_renv = [];
ripplesettings = [];
for i = 1:length(files)
    sessNum = num2str(files(i,1));
    trackInfo = num2str(files(i, 2));
    load([saveNeuralPath  '\' sessDate '_' sessNum '_'  trackInfo '\rawDataBySessionNeural.mat'])
    ripple_env = rawDataBySessionNeural.ripple.env;
    full_renv = [full_renv ripple_env];
end
for j = 1:size(full_renv, 1)
    ripplesettings(j).base = mean(full_renv(j,:));
    ripplesettings(j).stdev = std(full_renv(j,:));
end
for i = 1:length(files)
    sessNum = num2str(files(i,1));
    trackInfo = num2str(files(i, 2));
    savepath = fullfile(saveNeuralPath, [sessDate '_' sessNum '_'  trackInfo],'ripplesettings.mat');
    save(savepath, "ripplesettings")
end
end