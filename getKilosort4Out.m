function getKilosort4Out(subj, sessDate, neuralRawDataPath, dirs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% create single AP .bin file for all recordings on a given day %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Method 1:
% Steps:
%  1) Copy the .ap.bin file for each recording into a single location
%  2) Open command window and cd to that location
%  3) Use "COPY /B file1.bin + file2.bin file3.bin"
%    (e.g., COPY /B Neuropixels_practice_JK13_250425_2_g0_t0.imec0.ap.bin + Neuropixels_practice_JK13_250425_3_g0_t0.imec0.ap.bin +...
%    Neuropixels_practice_JK13_250425_4_g0_t0.imec0.ap.bin Neuropixels_practice_JK13_250425_merged_g0_t0.imec0.ap.bin)
%  https://stackoverflow.com/questions/53279744/how-to-join-two-binary-files-on-windows

%Method 2:
% Alternative, Automated Code to create *_merged_* files
tmpDirs = dir([neuralRawDataPath '\' sprintf('*%s_%s*',subj,sessDate)]);
tmp_ap_binName = sprintf([tmpDirs(1).name '_t0.imec0.ap.bin']);
fileName = [neuralRawDataPath '\' tmpDirs(1).name(1:end-4) 'merged_' tmp_ap_binName(end-17:end)];%session folder
altFileName = [neuralRawDataPath '\kilosort4\' tmpDirs(1).name(1:end-4) 'merged_' tmp_ap_binName(end-17:end)];%kilosort folder in session folder
outputFile = [];
if ~isfile(fileName) && ~isfile(altFileName)%if not a file
    for r = 1:size(tmpDirs,1)
        ap_binName = sprintf([tmpDirs(r).name '_t0.imec0.ap.bin']);
        ap_path = [tmpDirs(r).folder '\' tmpDirs(r).name '\' sprintf([tmpDirs(r).name '_imec0'])];
        if r == 1
            outputFile = fullfile(neuralRawDataPath, [tmpDirs(r).name(1:end-4) 'merged_' ap_binName(end-17:end)]);
            fid_out = fopen(outputFile, 'w');
        end
        fid_in = fopen([ap_path '\' ap_binName], 'r');
        while ~feof(fid_in)
            mergeData = fread(fid_in, 1e6, '*uint8');  % read in chunks
            fwrite(fid_out, mergeData, 'uint8');
        end
        fclose(fid_in);
    end
    fclose(fid_out);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% run kilosort4 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%generates kilosort4 folder with essential .npy files
%FOR PHY2: move merged ap .bin file into kilosort4 folder

%Method 1:
% Anaconda prompt:
%  1) Open GUI by typing 1. Y:; 2. activate kilosort; 3. cd 'folder with merged ap bin file'; 4. python -m kilosort
%  2) Run kilosort by clicking 1. "Select Binary File"; 2. "LOAD"; 3. "Run"
%  https://github.com/MouseLand/Kilosort

%Method 2:
% Auto-run kilosort4 in python. (if use kilosort3/2, use Nuri's master_kilosort2_defaultNJ.m)
if ~isfolder([tmpDirs(1).folder '\kilosort4'])
    pyenv('Version', dirs.kilosortPyEnv);
    if ~exist(fullfile(neuralRawDataPath, [tmpDirs(1).name(1:end-4) 'merged_' tmp_ap_binName(end-17:end)]), 'file')
        % not recommended. should create merged file before using kilosort.
        % post-processing steps in Kilosort4_Pipeline_XZ.m (or JLK version)
        % filter cells based on overall recording SNR etc
        sess_str = ['Neuropixels_' subj '_' sessDate];
        cmd = sprintf('"%s" %s --data_dir "%s" --sess_str "%s"', dirs.kilosortPyEnv, dirs.kilosortPyScript, neuralRawDataPath, sess_str);
    else
        % recommended method. use merged ap .bin files from the previous step
        cmd = sprintf('"%s" %s --data_dir "%s" --combined_ap_file "%s"', dirs.kilosortPyEnv, dirs.kilosortPyScript, neuralRawDataPath, [tmpDirs(1).name(1:end-4) 'merged_' tmp_ap_binName(end-17:end)]);
    end
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