function [files] = getfilenums(anrawdatadir)
%gets list of all available recording files for kilosort INTAN
%   need raw data folder
% Get all .rhd files in the directory
filestruct = dir(fullfile(anrawdatadir, 'recording*.rhd'));

% Extract all filenames
fnames = {filestruct.name};

% Preallocate
fileNums = zeros(1, numel(fnames));

% Loop through files and extract the number after "recording"
for i = 1:numel(fnames)
    % Example filename: recording1_251117_103439.rhd
    % Extract the digits after "recording" and before the first underscore
    tokens = regexp(fnames{i}, '^recording(\d+)_', 'tokens');
    fileNums(i) = str2double(tokens{1}{1});
end
files = unique(fileNums)
end