function [filteredData, y, time, timeNorm, timeDatetime] = preprocessData(ID, dataPos, data, time, fs_ori, fs_res, timeNorm, timeDatetime)
% preprocessData Applies filtering and resampling to input data.
%
% This function filters the input data with a low-pass filter and resamples
% it to the specified frequency. It also adjusts timestamps accordingly.
%
% INPUTS:
%   ID: Participant ID as a string (e.g., "23").
%   dataPos: Sensor position ('Body' or 'Wrist').
%   data: Raw input data (matrix or vector) to be filtered and resampled.
%   time: Timestamps corresponding to the input data.
%   timeNorm: Normalized time vector for the data.
%   timeDatetime: Datetime vector for the data.
%   fs_res: Target resampling frequency (Hz).
%
% OUTPUTS:
%   filteredData: Low-pass filtered input data.
%   y: Resampled version of the input data.
%   time: Resampled time vector.
%   timeNorm: Resampled normalized time vector.
%   timeDatetime: Resampled datetime vector.

% Define filter coefficients and target frequency
[LPa, LPb] = butter(4, 20/(fs_ori/2), 'low'); % Low-pass filter at 20 Hz

% Step 1: Apply the low-pass filter
filteredData = filtfilt(LPa, LPb, data);

% Step 2: Resample data and associated time vectors
% Resample data
[y, time] = resample(filteredData, time, fs_res);

% Resample normalized time
[~, timeNorm] = resample(filteredData, timeNorm, fs_res);
%timeNorm = time - time(1);

% Resample datetime vector
[~, timeDatetime] = resample(filteredData, timeDatetime, fs_res);

datause=y(:,4);
datausexyz=y(:,1:3);
% Save the processed data with aligned variable names
matName = strcat('ID', ID, '_processed.mat');
filename = strcat('dataProcessed_',dataPos);
if ~isfolder(filename)
    % If the folder does not exist, create it
    mkdir(filename);
    fprintf('Folder created: %s\n', filename);

end
save(fullfile(filename ,matName), 'datause','datausexyz', 'time', 'timeNorm', 'timeDatetime');
% Save the extracted features to the file
fprintf('Saved preprocessed data for ID%s in %s.\n' , ID, fullfile(fullfile(filename ,matName),matName));
end