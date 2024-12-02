%% Train models using HAR dataset - STEP 1: Read and preprocess data  
fs = 50;  % Sampling frequency in Hz
PeriodLimits = seconds([0.2, 7]);
WalkingPeriod = 9:18;  % CWT periods for walking activity
MaxCoefThres = 0.1;  % Threshold for CWT coefficients (max)
MaxCoefThres2 = 0.05;  % Threshold for CWT coefficients (mean)
MinWalkingPeriod = 5;  % Minimum walking period in seconds
ParID = [(1:18) + 500, 6, 8:10, 12:14, 16:20, 22:29];  % Participant IDs

% Paths
pathname70 = "Public Data\har70plus";
pathname = "Public Data\harth";

% Filter design
[LPa, LPb] = butter(4, 20 / (fs / 2), 'low');

% Preallocate data structures
filtMagdata = struct;
filtxyzdata = struct;
cwtCoef = struct;
Labels = struct;

NoParticipants = length(ParID);

% Loop through participants to read and preprocess data
for i = 1:NoParticipants
    % Select filename based on participant ID
    if i < 19
        filename = strcat(num2str(ParID(i)), '.csv');
        pathnamesel = pathname70;
    elseif ParID(i) < 10
        filename = strcat('S00', num2str(ParID(i)), '.csv');
        pathnamesel = pathname;
    else
        filename = strcat('S0', num2str(ParID(i)), '.csv');
        pathnamesel = pathname;
    end
    
    % Read data
    data = readmatrix(fullfile(pathnamesel, filename));
    
    % Filter magnitude and XYZ data
    filtdata = filtfilt(LPa, LPb, vecnorm(data(:, 2:4), 2, 2));  % Vector magnitude
    filtdataxyz = filtfilt(LPa, LPb, data(:, 2:4));  % Filtered XYZ components

    % Perform CWT with PeriodLimits and store results
    fprintf('Processing CWT for subject %d\n', i);
    [wceof, period] = cwt(filtdata - 1, seconds(1 / fs), 'PeriodLimits', PeriodLimits, 'amor');
    cwtCoef.(sprintf('Subject%d', i)) = wceof;
    filtMagdata.(sprintf('Subject%d', i)) = filtdata;
    filtxyzdata.(sprintf('Subject%d', i)) = filtdataxyz;
    Labels.(sprintf('Subject%d', i)) = data(:, end);
    fprintf('Finished subject %d\n', i);
end

%% STEP 2: Extract features
addpath('Functions')
window = 10;  % Window size in seconds
lyingOrient = 85;  % Lying orientation in degrees
staMethod=5; % Statistics method

for i = 1:NoParticipants
    fprintf('Extracting features for subject %d\n', i);
    
    wceof = cwtCoef.(sprintf('Subject%d', i));
    datause = filtMagdata.(sprintf('Subject%d', i));
    datausexyz = filtxyzdata.(sprintf('Subject%d', i));
    
    % Estimate body orientation
    BodyOrientEst=abs(rad2deg(atan(sqrt(datausexyz(:,2).^2+datausexyz(:,3).^2)./datausexyz(:,1))));
    ChangeOrient=abs(BodyOrientEst-circshift(BodyOrientEst,1));
    ChangeOrient(1)=0;
    
    % Find stationary periods using CWT coefficients
    stationary = and(max(abs(wceof(WalkingPeriod,:))) > MaxCoefThres,...
        mean(abs(wceof(WalkingPeriod,:))) > MaxCoefThres2);
    
    % Identify edges of stationary periods and remove short sequences
    d_data = diff([0 stationary 0]);
    startIdx = find(d_data == 1);  % Start indices of stationary periods
    endIdx = find(d_data == -1) - 1;  % End indices
    lengths = endIdx - startIdx + 1;
    
    % Remove short walking periods
    % shortSeqIdx = lengths < MinWalkingPeriod * fs;
    % stationary(startIdx(shortSeqIdx):endIdx(shortSeqIdx)) = 0;
    shortSeqIdx = find(lengths < MinWalkingPeriod*fs);

    % Convert the short 1 sequences to 0s
    for j = 1:length(shortSeqIdx)
        stationary(startIdx(shortSeqIdx(j)):endIdx(shortSeqIdx(j))) = 0;
    end
    startIdx(shortSeqIdx)=[];
    endIdx(shortSeqIdx)=[];

    % Define windowed feature extraction
    idx = 1:(window * fs):length(datause);
    if idx(end) ~= length(datause)
        idx = [idx, length(datause)];
    end
    
    % Extract features from each window
    for j = 1:(length(idx) - 1)
        timerange = idx(j):idx(j+1);
        if j==1
            ChangeOrient1=ChangeOrient(timerange(2:end));
        else
            ChangeOrient1=ChangeOrient(timerange);
        end

        BodyOrientEstuse = BodyOrientEst(timerange);

        absCoef = abs(wceof(:, timerange));
        absCoefWalking = abs(wceof(WalkingPeriod, timerange));

        % Feature extraction using peaks in the walking coefficient
        [peaks1, peaks1Loc] = findpeaks(mean(absCoefWalking), 'MinPeakHeight', 0.05, 'MinPeakProminence', 0.01);

        if isempty(peaks1)
            staout1 = [zeros(1,8),0];
        else
            staout1 = [staset(peaks1', 1, 6), staset(peaks1Loc', 1, 6), length(peaks1)];
        end

        [maxv_pr, maxvl_pr] = max(max(absCoef, [], 2));
        [maxv_wr, maxvl_wr] = max(max(absCoefWalking, [], 2));

        % Aggregate features
        staout = [
            staset(datause(timerange), 1, staMethod), ...
            staset(datausexyz(timerange, 1), 1, staMethod), ...
            staset(max(absCoef)', 1, staMethod), ...
            staset(mean(absCoef)', 1, staMethod), ...
            staset(max(absCoefWalking)', 1, staMethod), ...
            staset(mean(absCoefWalking)', 1, staMethod), ...
            staset(mean(absCoefWalking)' ./ mean(absCoef)', 1, staMethod), ...
            staout1, ...
            time2num(period(maxvl_pr)), maxv_pr, time2num(period(maxvl_wr)), maxv_wr, ...
            staset(ChangeOrient1, 1, staMethod), ...
            staset(BodyOrientEstuse, 1, staMethod), ...
            max(BodyOrientEstuse) - min(BodyOrientEstuse), ...
            mean(BodyOrientEstuse) - lyingOrient, ...
            sum(stationary(timerange))
        ];

        % Store labels and features
        labelint = mode(Labels.(sprintf('Subject%d', i))(timerange), 'all');
        wlabelint = mode(stationary(timerange), 'all');
        
        if and(i==1,j==1)
            RealLabel=labelint;
            RealWalkingLabel=(labelint==1||labelint==4||labelint==5);
            PredictWalkingLabel=wlabelint;
            ParticipantLabel= i;
            f=staout;
        else
            RealLabel=[RealLabel;labelint];
            RealWalkingLabel=[RealWalkingLabel;labelint==1];
            ParticipantLabel=[ParticipantLabel;i];
            PredictWalkingLabel=[PredictWalkingLabel;wlabelint];
            f=[f;staout];
        end
    end
    fprintf('Finished extracting features for subject %d\n', i);
end
%% Define labels
% standing with leg movements as standing
WalkingSSLO=string(size(RealLabel));
WalkingSSLO(any([(RealLabel==1),(RealLabel==4),(RealLabel==5)],2))="Walking";
WalkingSSLO(any([(RealLabel==3),(RealLabel==6)],2))="Standing";
WalkingSSLO(RealLabel==7)="Sitting";
WalkingSSLO(RealLabel==8)="Lying";
% WalkingSSLO(RealLabel==2)="Running";
WalkingSSLO(any([(RealLabel==2),(RealLabel==13),(RealLabel==14),(RealLabel==130),(RealLabel==140)],2))="Other";
%% Step3: save to csv files for model training
% Define feature names
feature_names = ["AccMag_mean","AccMag_std","AccMag_var","AccMag_median",...
                  "AccMag_25perc","AccMag_75perc","StandAxis_mean","StandAxis_std",...
                  "StandAxis_var","StandAxis_median","StandAxis_25perc",...
                  "StandAxis_75perc","ScMaxPAllf_mean","ScMaxPAllf_std","ScMaxPAllf_var",...
                  "ScMaxPAllf_median","ScMaxPAllf_25perc","ScMaxPAllf_75perc",...
                  "ScMeanPAllf_mean","ScMeanPAllf_std","ScMeanPAllf_var","ScMeanPAllf_median",...
                  "ScMeanPAllf_25perc","ScMeanPAllf_75perc","ScMaxPWalkf_mean",...
                  "ScMaxPWalkf_std","ScMaxPWalkf_var","ScMaxPWalkf_median","ScMaxPWalkf_25perc",...
                  "ScMaxPWalkf_75perc","ScMeanPWalkf_mean","ScMeanPWalkf_std","ScMeanPWalkf_var",...
                  "ScMeanPWalkf_median","ScMeanPWalkf_25perc","ScMeanPWalkf_75perc",...
                  "ScRatioP_mean","ScRatioP_std","ScRatioP_var","ScRatioP_median","ScRatioP_25perc",...
                  "ScRatioP_75perc","CoefPeaks_mean","CoefPeaks_std","CoefPeaks_var",...
                  "CoefPeaks_median","CoefPeakWidths_mean","CoefPeakWidths_std","CoefPeakWidths_var",...
                  "CoefPeakWidths_median","CoefPeaks_","DomfAll_","TiAveDomfAll_","DomfWalk_","TiAveDomfWalk_",...
                  "OrientChange_mean","OrientChange_std","OrientChange_var","OrientChange_median",...
                  "OrientChange_25perc","OrientChange_75perc","Orient_mean","Orient_std","Orient_var",...
                  "Orient_median","Orient_25perc","Orient_75perc","MaxOrientChange_","MeanOrientDiffLying_","Perctime_"];

% Create a table for the features
features_table = array2table(f, 'VariableNames', cellstr(feature_names));

% Add labels and participant information
features_table = [features_table,table(WalkingSSLO', 'VariableNames', {'Labels'})];

% Update WalkingSSLO for combined 'Standing' and 'Sitting' labels
WalkingSSLO(ismember(WalkingSSLO, {'Standing', 'Sitting'})) = {'Standing/Sitting'};
T=[table(WalkingSSLO','VariableNames', {'Labels2'}),table(ParticipantLabel,'VariableNames', {'ParLabel'})];
features_table = [features_table,T];

% Save the table to a CSV file
if ~exist('features', 'dir')
    mkdir('features');  % Create the folder if it doesn't exist
end
writetable(features_table, fullfile('features', 'processed_features_matlab.csv'));

disp("Features saved to 'features\processed_features.csv'");
%% Helper function
function [staout]=staset(V,dim,method)
switch method
    case 1
        staout=[mean(V,dim),std(V,dim),(std(V,dim)./mean(V,dim))*100,median(V,dim),...
            skewness(V,dim),kurtosis(V,dim)];
    case 2
        staout=[mean(V,dim),std(V,dim),(std(V,dim)./mean(V,dim))*100,median(V,dim),...
          max(V,[],dim),min(V,[],dim)];

    case 5
        staout=[mean(V,dim),std(V,dim),(std(V,dim)./mean(V,dim))*100,median(V,dim),...
            prctile(V,25,dim),prctile(V,75,dim)];
    case 6
        staout=[mean(V,dim),std(V,dim),(std(V,dim)./mean(V,dim))*100,median(V,dim)];
end
end