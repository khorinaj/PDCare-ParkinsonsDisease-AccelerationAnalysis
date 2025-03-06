%% addpath
addpath("utils\")

%% ****************************************************************************************************
% *****************************************************************************************************
%******************************************************************************************************

%                 Process a set of PD data using TRUNK sensor


%******************************************************************************************************
% *****************************************************************************************************
% *****************************************************************************************************
IDs=["23","27","05","16","07","09","06","01"];
dataPosition="Body"; % wrist-2, body-1
%%
fs_ori = 100; % original sampling rate
fs_res = 50; % Target resampling frequency (Hz)

% Specify the path of raw data
dataDir = 'C:\Users\zceeysu\OneDrive - University College London\PD long term NEW';
windowSize=10; %s windon size for feature extractions

for i=1:length(IDs)

    ID=IDs(i);
    % % Preprocessing PD data
    % %Load the participant data
    % load(fullfile(dataDir,strcat('ID',ID)), 'dataBodyINT', 'time', 'timeNorm', 'timeDatetime');
    % % Preprocess the data
    % preprocessData(ID,dataPosition,dataBodyINT, time, fs_ori,  fs_res, timeNorm, timeDatetime);

    % % Extract features
    % % define axis pointing to gravity when standing
    % load('Info\standAxis','standAxis') % axis index
    % load('Info\gdirection','gdirection') % axis direction
    % fPD = extractFeaturesBody(ID,standAxis.(strcat('ID',ID)), gdirection.(strcat('ID',ID)),windowSize);

    % Extract labels predicted by python model
    % filename for storing the predicted labels
    filename=fullfile(strcat('dataProcessed_',dataPosition),strcat('ID',ID,'Label4classfs_pt2.xlsx'));
    label=extractLabelsBody(ID,windowSize,filename);
    % 
    % Visualize activity predicted
    fontsize=30;
    activityVisualizeBody(ID,fontsize)
end

%% Extract activity features based on labels predicted by python model
IDstrings=IDs;
Chunkslength="24h";
WalkingLengthLim=10;
OT=5;
[T,fcorrALL]=...
    extractActivityFeatures(IDstrings,dataPosition,Chunkslength,WalkingLengthLim,'OutlierThreshold',OT);

%% Correlation Analysis
methodsel="Spearman";
methodsel="Kendall";
T=T(table2array(T(:,2))>21.6,:);
% SL= PDadditional3([3,6,7,4,1,2] , :);
% SL=readtable('clinicalScales\ClinicalScales.xlsx','Sheet','UPDRS_TRUNK');
% scaleName='UPDRS';
% SL=readtable('clinicalScales\ClinicalScales.xlsx','Sheet','Other_TRUNK');
% scaleName='Total';
SL=readtable('clinicalScales\ClinicalScales.xlsx','Sheet','Paper_TRUNK');
scaleName='Paper';
[values1,pvtable,Tfeature]=correlationAnalysis(T,dataPosition,methodsel,SL,scaleName, false);


%% ****************************************************************************************************
% *****************************************************************************************************
%******************************************************************************************************

%                 Process a set of PD data using WRIST sensor


%******************************************************************************************************
% *****************************************************************************************************
% *****************************************************************************************************

IDs=["23","27","05","16","19","18"];
dataPosition="Wrist"; % wrist-2, body-1
%%
fs_ori = 100; % original sampling rate
fs_res = 50; % Target resampling frequency (Hz)

% Specify the path of raw data
dataDir = 'C:\Users\zceeysu\OneDrive - University College London\PD long term NEW';
windowSize=10; %s windon size for feature extractions
for i=1:1
%for i=1:length(IDs)
    ID=IDs(i);
    % Preprocessing PD data
    % Load the participant data
    load(fullfile(dataDir,strcat('ID',ID)), 'dataWristINT', 'time', 'timeNorm', 'timeDatetime');
    % Preprocess the data
    preprocessData(ID,dataPosition, dataWristINT, time, fs_ori,  fs_res, timeNorm, timeDatetime);
end
%% extract labels predicted by biobank models
dataDir='C:\Users\zceeysu\OneDrive - University College London\PD long term NEW\GitHub Paper\PDdata\csvWristData';

ClassificationType="ActivityType"; % "ActivityType" or "ActivityIntensity";
extractLabelBiobank(IDs,ClassificationType,dataDir)

ClassificationType="ActivityIntensity"; % "ActivityType" or "ActivityIntensity";
extractLabelBiobank(IDs,ClassificationType,dataDir)

%% Extract activity features based on labels predicted by Biobank model
IDstrings=IDs;
Chunkslength="24h";
WalkingLengthLim=10;
OT=5;
[T_wrist,fcorrALL]=...
    extractActivityFeatures(IDstrings,dataPosition,Chunkslength,WalkingLengthLim,'OutlierThreshold',OT);

%% Correlation Analysis
methodsel="Spearman";
%methodsel="Kendall";
T=T(table2array(T(:,2))>21.6,:);
% SL=readtable('clinicalScales\ClinicalScales.xlsx','Sheet','UPDRS_WRIST');
% scaleName='UPDRS';
% SL=readtable('clinicalScales\ClinicalScales.xlsx','Sheet','Other_WRIST');
% scaleName='Total';
SL=readtable('clinicalScales\ClinicalScales.xlsx','Sheet','Paper_WRIST');
scaleName='Paper';
[values1_wrist,pvtable_wrist,Tfeature_wrist]=correlationAnalysis(T_wrist,dataPosition,methodsel,SL,scaleName, false);
