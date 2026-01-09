function getKilosort4Out_intan(subj, sessDate, neuralRawDataPath, dirs, files, params)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% create single .bin file for all recordings on a given day %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%making .bin file based on Kilosort2_Pipeline, rest based on
%getKilosort4Out from JK

%input animal info and project settings
addpath('Y:\singer\Danielle\code\kilosort2-pipeline\kilosort2-pipeline\functions');
%[params, dirs, th] = userProfiles_K2pipeline('Danielle', 'DAlesion');%go back and lump in correct from the main pipeline
%%%update everything here to match other formatting%%%

%get the animal info based on the inputs
addpath('Y:\singer\Danielle\code\preprocessing_Xv\preprocessing-pipeline\Functions');
sessions = sessDate; %define session as one date
%NN probe channel map - Takahashi_Intan_um_kilosortChanMap
%CN probe channel map - CN_doublesided_P2D_KSchmap
% anrawdatadir = neuralRawDataPath;
run.preCuration = 0;            %write specificed files to .bin for Kilosort

%% Run Kilosort Pipeline
for d = 1:size(sessions,1)
    %% get the animal/day info
    %sessionInfo = allindex(ismember(allindex(:,1:2), sessions(d,:), 'rows'),:);
    params.files = files;
    params.animal = subj;
    params.day =  sessDate;
    params.brainReg = 'CA1';%allindexT{allindexT.Animal == params.animal & allindexT.Date == params.day & allindexT.Recording == params.files(1),{'RegAB','RegCD'}};
    anrawdatadir = [dirs.rawdata, params.animal, '_', num2str(params.day), '\'];
    tempfiledir = [dirs.processeddata, params.animal, '_', num2str(params.day), '\'];
    anclusterdir = [dirs.cluster_local, params.animal, '_', num2str(params.day), '\'];
    %% write raw recording files to BIN for kilosort2
    if run.preCuration
        if ~exist(anclusterdir, 'dir'); mkdir(anclusterdir); end
        converttoBIN_K2(anrawdatadir, anclusterdir, params.files, params.probeChannels, params.brainReg, dirs.clusfolder)
    end
end

%%now bin made, actually run kilosort

%Method 2:
% Auto-run kilosort4 in python. (if use kilosort3/2, use Nuri's master_kilosort2_defaultNJ.m)
kilosortFolder = fullfile(anclusterdir, 'CA1', 'sorted', 'kilosort');
if ~isfolder([kilosortFolder '\kilosort4'])
    pyenv('Version', dirs.kilosortPyEnv);
    combinedFile = fullfile(kilosortFolder, 'allrecordings.bin');
    cmd = sprintf('"%s" %s --data_dir "%s" --combined_ap_file "%s"', dirs.kilosortPyEnv, dirs.kilosortPyScript, kilosortFolder, combinedFile);
    tic
    disp('Running kilosort4...')
    [status, result] = system(cmd);
    toc
    disp('Finish running kilosort4')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% manually spike sort the kilosort output using phy2 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Method 1:
% Anaconda prompt:
%  1) Open GUI by typing 1. Y:; 2. activate phy2; 3. cd 'folder with merged ap bin file'; 4. phy template-gui params.py

end%function