function fPD = extractFeaturesBody(ID,standAxis,gdirection,window)
% extractFeaturesForID Extract features for a single participant ID.
%
% This function loads preprocessed data for a specific participant,
% computes wavelet and statistical features over specified windows, and
% returns the extracted features.
%
% INPUTS:
%   ID: Participant ID as a string (e.g., "23").
%   standAxis: Index for the standing axis (pointing to gravity when standing)
%              used for body orientation.
%   gdirection: Direction multiplier for the orientation.
%   window: Window length in seconds for feature extraction.
%
% OUTPUTS:
%   fPD: Extracted features for the participant.


% Parameters for feature extraction
PeriodLimits = seconds([0.2, 7]); % CWT period limits
WalkingPeriod = 9:18; % Indices for walking periods (0.35s-0.65s)
MaxCoefThres = 0.1; % Threshold for maximum wavelet coefficients
MaxCoefThres2 = 0.05; % Threshold for mean wavelet coefficients
staMethod = 5; % Statistical aggregation method
lyingOrient = 85; % Lying orientation threshold
MinWalkingPeriod = 5; % Minimum walking period in seconds

% Load preprocessed data
matName = strcat('ID', ID, '_processed');
load(fullfile('dataProcessed_Body',matName), 'datause', 'datausexyz', 'timeNorm');
fs = 1 / (timeNorm(2) - timeNorm(1)); % Calculate sampling frequency (Hz)

% Calculate standing axis
StandAxis = gdirection * datausexyz(:, standAxis);

% Perform Continuous Wavelet Transform (CWT) on the data
[wceof, period, ~] = cwt(datause - 1, seconds(1 / fs), 'PeriodLimits', PeriodLimits, 'amor');

% Calculate body orientation estimate based on accelerometer data
BodyOrientEst = abs(rad2deg(atan(sqrt(datausexyz(:, 2).^2 + datausexyz(:, 3).^2) ./ datausexyz(:, 1) * gdirection)));

% Handle special cases for participant IDs
if ID == "27"
    load('Info/DirChange', 'DirChange');
    DirChange = DirChange.(strcat('ID', ID));
    StandAxis(DirChange.DataIndex:end) = datausexyz(DirChange.DataIndex:end, 2) * gdirection;
    BodyOrientEst(DirChange.DataIndex:end) = abs(rad2deg(atan(sqrt(datausexyz(DirChange.DataIndex:end, 1).^2 + datausexyz(DirChange.DataIndex:end, 3).^2) ./ datausexyz(DirChange.DataIndex:end, 2) * gdirection)));
elseif ID == "16"
    load('Info/DirChange', 'DirChange', 'DirChange2');
    ChangeDirLoc = DirChange.(strcat('ID', ID));
    ChangeDirLocS = DirChange2.(strcat('ID', ID));
    StandAxis(ChangeDirLocS.DataIndex:ChangeDirLoc.DataIndex) = datausexyz(ChangeDirLocS.DataIndex:ChangeDirLoc.DataIndex, 3);
    BodyOrientEst(ChangeDirLocS.DataIndex:ChangeDirLoc.DataIndex) = abs(rad2deg(atan(sqrt(datausexyz(ChangeDirLocS.DataIndex:ChangeDirLoc.DataIndex, 1).^2 + datausexyz(ChangeDirLocS.DataIndex:ChangeDirLoc.DataIndex, 2).^2) ./ datausexyz(ChangeDirLocS.DataIndex:ChangeDirLoc.DataIndex, 3))));
end

% Compute the change in orientation
ChangeOrient = abs(BodyOrientEst - circshift(BodyOrientEst, 1));
ChangeOrient(1) = 0; % Correct the first element

% Detect stationary periods based on wavelet coefficients
stationary = and(max(abs(wceof(WalkingPeriod, :))) > MaxCoefThres, ...
    mean(abs(wceof(WalkingPeriod, :))) > MaxCoefThres2);

% Find the start and end indices of stationary sequences
d_data = diff([0 stationary 0]);
startIdx = find(d_data == 1); % Start indices of stationary sequences
endIdx = find(d_data == -1) - 1; % End indices of stationary sequences

% Filter out short stationary sequences
lengths = endIdx - startIdx + 1;
shortSeqIdx = find(lengths < MinWalkingPeriod * fs);

% Remove short sequences from the stationary signal
for j = 1:length(shortSeqIdx)
    stationary(startIdx(shortSeqIdx(j)):endIdx(shortSeqIdx(j))) = 0;
end

% Recompute the start and end indices of stationary sequences
startIdx(shortSeqIdx) = [];
endIdx(shortSeqIdx) = [];

% Define processing windows
idx = 1:(window * fs):length(datause);
if idx(end) ~= length(datause)
    idx = [idx, length(datause)];
end

% Initialize the feature matrix
fPD = zeros((length(idx) - 1), 70);

% Extract features for each window
for j = 1:(length(idx) - 1)
    timerange = idx(j):idx(j + 1);

    % Compute the change in orientation for the current window
    if j == 1
        ChangeOrient1 = ChangeOrient(timerange(2:end));
    else
        ChangeOrient1 = ChangeOrient(timerange);
    end

    % Extract the body orientation estimate for the current window
    BodyOrientEstuse = BodyOrientEst(timerange);

    % Compute wavelet coefficients for the current window
    absCoef = abs(wceof(:, timerange));
    absCoefWalking = abs(wceof(WalkingPeriod, timerange));

    % Perform wavelet peak detection
    [peaks1, peaks1Loc] = findpeaks(mean(absCoefWalking), 'MinPeakHeight', 0.05, 'MinPeakProminence', 0.01);

    % If no peaks are detected, assign default values
    if isempty(peaks1)
        staout1 = [zeros(1, 8), 0];
    else
        staout1 = [staset(peaks1', 1, 6), staset(peaks1Loc', 1, 6), length(peaks1)];
    end

    % Extract statistical features for the current window
    [maxv_pr, maxvl_pr] = max(max(absCoef, [], 2));
    [maxv_wr, maxvl_wr] = max(max(absCoefWalking, [], 2));

    staout = [staset(datause(timerange), 1, staMethod), ...
        staset(StandAxis(timerange), 1, staMethod), ...
        staset(max(absCoef)', 1, staMethod), ...
        staset(mean(absCoef)', 1, staMethod), ...
        staset(max(absCoefWalking)', 1, staMethod), ...
        staset(mean(absCoefWalking)', 1, staMethod), ...
        staset(mean(absCoefWalking)' ./ mean(absCoef)', 1, staMethod), ...
        staout1, ...
        time2num(period(maxvl_pr)), maxv_pr, ...
        time2num(period(maxvl_wr)), maxv_wr, ...
        staset(ChangeOrient1, 1, staMethod), ...
        staset(BodyOrientEstuse, 1, staMethod), ...
        max(BodyOrientEstuse) - min(BodyOrientEstuse), ...
        mean(BodyOrientEstuse) - lyingOrient, ...
        sum(stationary(timerange))];

    % Store features in the feature matrix
    fPD(j, :) = staout;
end
% feature_names = ["AccMag_mean","AccMag_std","AccMag_var","AccMag_median",...
%                   "AccMag_25perc","AccMag_75perc","StandAxis_mean","StandAxis_std",...
%                   "StandAxis_var","StandAxis_median","StandAxis_25perc",...
%                   "StandAxis_75perc","ScMaxPAllf_mean","ScMaxPAllf_std","ScMaxPAllf_var",...
%                   "ScMaxPAllf_median","ScMaxPAllf_25perc","ScMaxPAllf_75perc",...
%                   "ScMeanPAllf_mean","ScMeanPAllf_std","ScMeanPAllf_var","ScMeanPAllf_median",...
%                   "ScMeanPAllf_25perc","ScMeanPAllf_75perc","ScMaxPWalkf_mean",...
%                   "ScMaxPWalkf_std","ScMaxPWalkf_var","ScMaxPWalkf_median","ScMaxPWalkf_25perc",...
%                   "ScMaxPWalkf_75perc","ScMeanPWalkf_mean","ScMeanPWalkf_std","ScMeanPWalkf_var",...
%                   "ScMeanPWalkf_median","ScMeanPWalkf_25perc","ScMeanPWalkf_75perc",...
%                   "ScRatioP_mean","ScRatioP_std","ScRatioP_var","ScRatioP_median","ScRatioP_25perc",...
%                   "ScRatioP_75perc","CoefPeaks_mean","CoefPeaks_std","CoefPeaks_var",...
%                   "CoefPeaks_median","CoefPeakWidths_mean","CoefPeakWidths_std","CoefPeakWidths_var",...
%                   "CoefPeakWidths_median","CoefPeaks_","DomfAll_","TiAveDomfAll_","DomfWalk_","TiAveDomfWalk_",...
%                   "OrientChange_mean","OrientChange_std","OrientChange_var","OrientChange_median",...
%                   "OrientChange_25perc","OrientChange_75perc","Orient_mean","Orient_std","Orient_var",...
%                   "Orient_median","Orient_25perc","Orient_75perc","MaxOrientChange_","MeanOrientDiffLying_","Perctime_"];
% 
% % Create a table for the features
% fPD = array2table(fPD, 'VariableNames', cellstr(feature_names));
save(fullfile('dataProcessed_Body',matName), "fPD", '-append');
% Save the extracted features to the file
fprintf('Saved feature extraction for ID%s in %s.\n' , ID, fullfile('dataProcessed_Body',matName));
end
