function  [T,fcorrALL]=extractActivityFeatures(ParID,dataPosition,Chunkslength,WalkingLengthLim,varargin)
% extractActivityFeatures - Extract activity features from long-term sensor data
%
% Description:
% This function processes long-term data collected from trunk or wrist sensors and 
% extracts activity features based on predicted activity labels. The data is split 
% into chunks (e.g., 24h, 48h) for each participants to allow analysis over specific time intervals.
% It provides both a detailed feature table (`T`) and a matrix of extracted features (`fcorrALL`).
%
% Inputs:
% - ParID: (array of strings or numbers) Participant IDs to process (e.g., ["01", "02"]).
% - dataPosition: (string) Specifies the sensor position:
%     - 'Body': Data from trunk sensors.
%     - 'Wrist': Data from wrist sensors.
% - Chunkslength: (string) Defines the chunk size for analysis:
%     - '24h': Splits data into 24-hour intervals.
%     - '48h': Splits data into 48-hour intervals.
%     - '24h-fix': Splits data based on calendar days starting at 12:00 PM.
%     - 'all': Uses the entire dataset as a single chunk.
% - WalkingLengthLim: (numeric) Minimum threshold for walking length to consider.
% - varargin: (optional) Customizable parameters for the function:
%     - 'removeOutlierFlag': (boolean, default = true) Flag to remove outliers from walking peaks.
%     - 'OutlierThreshold': (numeric, default = 15) Z-score threshold for outlier detection.
%     - 'MinPeakDistance': (numeric, default = 0.35) Minimum peak distance for walking feature detection.
%
% Outputs:
% - T: (table) A table containing extracted features for each participant and chunk. Columns include:
%     - 'DataID': Participant IDs.
%     - 'NoDataSample(h)': Duration of data available in hours.
%     - Activity feature columns (specific to sensor position, e.g., mean walking speed).
% - fcorrALL: (matrix) A matrix of all extracted features, where each row corresponds to a participant-chunk combination.

% Key Features:
% - Processes data from trunk or wrist sensors, adapting to the sensor type.
% - Splits data into customizable chunks for time-based analysis (e.g., daily or entire dataset).
% - Extracts a comprehensive set of activity features, including:
%     - Walking-related features (e.g., peak walking speed, duration, gaps between walking bouts).
%     - Statistical summaries of activity data (mean, standard deviation, etc.).
% - Handles outliers in walking-related features with customizable thresholds.
% - Outputs features in a structured table format for easy downstream analysis.

% Example Usage:
% [T, fcorrALL] = extractActivityFeatures(["01", "02"], 'Body', '24h', 0.5, 'removeOutlierFlag', true, 'OutlierThreshold', 10);
% This extracts activity features from trunk sensor data for participants '01' and '02', 
% splitting the data into 24-hour chunks and removing outliers with a Z-score threshold of 10.

% Dependencies:
% - Requires helper functions such as `findWalkFeature2`, `getBodyStatistics`, and `getWristStatistics`.
% - Loads sensor-specific activity feature names from "Info/ActivityFeatureName.mat".
% - Data files must be stored in `dataProcessed_Body/` or `dataProcessed_Wrist/`.



addpath("utils\wavelet-coherence-master\")

removeOutlierFlag=true;
for i = 1:numel(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i}, 'removeOutlierFlag')
        % Look for 'method' in varargin and set the method variable
        if i < numel(varargin)
            removeOutlierFlag = varargin{i + 1};
            break
        end

    end
end

OutlierThreshold=15;
for i = 1:numel(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i}, 'OutlierThreshold')
        % Look for 'method' in varargin and set the method variable
        if i < numel(varargin)
            OutlierThreshold = varargin{i + 1};
            break
        end

    end
end


MinPeakDistance=0.35;
for i = 1:numel(varargin)
    if ischar(varargin{i}) && strcmpi(varargin{i}, 'MinPeakDistance')
        % Look for 'method' in varargin and set the method variable
        if i < numel(varargin)
            MinPeakDistance = varargin{i + 1};
            break
        end

    end
end


DataID="1";
DataLength=0;
count=0;
dataDir=strcat('dataProcessed_',dataPosition);

for io=1:(length(ParID))
    %----------Load data and preparation --------------------------------------
    count=count+1;
    ID=(ParID(io));
    
    matName = strcat('ID', ID, '_processed');
    if dataPosition=="Body"
        load(fullfile(dataDir,matName), 'datause','timeDatetime','time','timeNorm','label')

    elseif dataPosition=="Wrist"
         load(fullfile(dataDir,matName),'datause','timeDatetime','time','timeNorm')
         load(fullfile(dataDir,strcat('ID',ID,'ATlabel')),'expandedLabels')
         label=expandedLabels(:,1:6);
         load(fullfile(dataDir,strcat('ID',ID,'AIlabel')),'expandedLabels')
         label=[label,expandedLabels(:,1:4)];
         label(end,:)=label(end-1,:);

    else
        error('Wrong sensor position')
    end


    fs=round(1/(timeNorm(2)-timeNorm(1)));
    Nodata24h=24*60*60*fs;
    %---------Split data to chunk--------------------------------------

    if Chunkslength=="48h"
        remainder = mod(length(datause), (Nodata24h*2));
        if remainder == 0
            NoChunks=length(datause)/(Nodata24h*2);
        else
            NoChunks=floor(length(datause)/(Nodata24h*2))+1;
        end

        no24h=2;
    elseif Chunkslength=="24h"
        remainder = mod(length(datause), (Nodata24h));

        if remainder == 0
            NoChunks=length(datause)/(Nodata24h);
        else
            NoChunks=floor(length(datause)/(Nodata24h))+1;
        end
        no24h=1;
    elseif Chunkslength=="24h-fix"
        % Find unique days in the timeDatetime
        days = unique(dateshift(time, 'start', 'day'));
        days=days+hours(12);
        NoChunks=(length(days)-1);

    elseif Chunkslength=="all"
        NoChunks=1;

    elseif Chunkslength=="all-fix"
        NoChunks=1;
        timenow=time(1);
        while timenow<time(end)
            timenow=timenow+day(1);
        end
        [~,endpoint]=min(abs(time-timenow));
    else
        error('Wrong Chunkslength input')
    end

    %---------------------------------------------------------------------------

    if dataPosition=="Body"
        fcorr=zeros(NoChunks,23);

    elseif dataPosition=="Wrist"
        fcorr=zeros(NoChunks,29);
    else
        error('Wrong sensor position')
    end
    % fcorr=zeros(NoChunks,23);
    tstartIdx=1;
    for ii=1:NoChunks

        if Chunkslength=="24h-fix"
            RangeData=(time >= days(ii)) & (time <= days(ii+1));
            datalen=sum(RangeData);
        else

            if or(Chunkslength=="48h",Chunkslength=="24h")
                if tstartIdx+(Nodata24h*no24h)<=length(datause)
                    tendIdx=tstartIdx+(Nodata24h*no24h);
                else
                    tendIdx=length(datause);
                end
            elseif Chunkslength=="all"
                tendIdx=length(datause);
            else
                tendIdx=endpoint;
            end
            RangeData=tstartIdx:tendIdx;
            datalen=length(tstartIdx:tendIdx);

            tstartIdx=tendIdx+1;

        end

        label24h=label(RangeData,:);
        normt24h=time(RangeData);
        time24h=timeDatetime(RangeData);
        datause24h=datause(RangeData);

        DataLength=[DataLength;datalen];
        DataID=[DataID;strcat('ID',ID)];
    
        if dataPosition=="Body"
            WalkingLabel=(label24h==1);

        elseif dataPosition=="Wrist"
            WalkingLabel=boolean(label24h(:,6));
        end


        if sum(WalkingLabel)==0
            Walkingfeatures=[nan(1,9)];
        else
            [WalkPeaks,WalkGaps,~,~]=findWalkFeature2(WalkingLabel,datause24h,normt24h,fs,MinPeakDistance,WalkingLengthLim,0);

            if removeOutlierFlag
                z_scores = zscore(WalkPeaks);
                % Identify outliers
                outliers = abs(z_scores) > 15;
                WalkPeaks=WalkPeaks(~outliers);
            end
            Walkingfeatures=[length(WalkPeaks),staset(WalkPeaks,1,6),staset(WalkGaps(WalkGaps<1.5),1,6)];
        end

        if dataPosition=="Body"
            fcorr(ii,:)=[getBodyStatistics(label24h,time24h,datause24h,removeOutlierFlag,OutlierThreshold),Walkingfeatures];

        elseif dataPosition=="Wrist"
            fcorr(ii,:)=[getWristStatistics(label24h,time24h,datause24h,removeOutlierFlag,OutlierThreshold),Walkingfeatures];

        end

    end

    sprintf('Finished: ID%s',ID)
    if io==1
        fcorrALL=fcorr;

    else
        fcorrALL=[fcorrALL;fcorr];
    end

end
disp(strcat(dataPosition ,'features extracted'))
ParID

DataID(1)=[];
DataLength(1)=[];


if dataPosition=="Body"
    load("Info\ActivityFeatureName.mat",'BodyfeatureName')
    featureName=BodyfeatureName;
elseif dataPosition=="Wrist"
    load("Info\ActivityFeatureName.mat",'WristfeatureName')
    featureName=WristfeatureName;
end
% Create the main table for fcorrALL
T = array2table([DataLength/Nodata24h*24, fcorrALL], 'VariableNames', [{'NoDataSample(h)'}, cellstr(featureName)]);

% Add DataID as the first column
T = [table(DataID, 'VariableNames', {'DataID'}), T];

end

