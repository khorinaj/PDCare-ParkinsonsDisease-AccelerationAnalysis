"""
Script Name: feature_extraction_pipeline.py

Description:
This script preprocesses the data in public dataset, performs feature extraction, and saves the extracted features
for machine learning model training. The key steps include:
1. Preprocessing: Filter raw sensor data, compute magnitude and XYZ components, and apply
   continuous wavelet transform (CWT) to extract coefficients.
2. Feature Extraction: Extract time- and frequency-domain features from the processed data
   using sliding windows and compute statistical features.
3. Labeling and Export: Map activity labels to predefined categories, and save the extracted
   features and labels to a CSV file for further analysis.

Inputs:
- Raw sensor data in CSV format stored in the "Public Data" directories.
- Participant IDs (ParID) to select appropriate datasets.

Outputs:
- Processed feature dataset saved as 'features/processed_features_python.csv'.

Usage:
Run the script as follows:
    python feature_extraction_pipeline.py
Ensure all dependencies are installed and input files are in the correct directories.

Author: Yanke Sun
Date: 31/10/2024
"""
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt,find_peaks
from scipy.stats import skew, kurtosis, mode
import pywt  # for wavelet transform (CWT)
import joblib

# Paths (change the file if available)
pathname70 = "Public Data/har70plus"
pathname = "Public Data/harth"
# Participant ID to precess
ParID = np.concatenate([np.arange(1, 19) + 500, [6, 8, 9, 10, 12, 13, 14, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29]])

# Define parameter for preprocessing
fs = 50  # Sampling frequency
# Butterworth filter design
LPa, LPb = butter(4, 20 / (fs / 2), 'low')

# Define parameter for cwt
wavelet='morl'
# cwt fequency limits 
frequencies=[5,4.66516495768404,4.35275281648062,4.06126198178118,3.78929141627600,3.53553390593274,3.29876977693224,3.07786103336229,2.87174588749259,2.67943365634073,2.50000000000000,2.33258247884202,2.17637640824031,2.03063099089059,1.89464570813800,1.76776695296637,1.64938488846612,1.53893051668115,1.43587294374629,1.33971682817037,1.25000000000000,1.16629123942101,1.08818820412016,1.01531549544530,0.947322854068999,0.883883476483185,0.824692444233060,0.769465258340573,0.717936471873148,0.669858414085184,0.625000000000001,0.583145619710505,0.544094102060078,0.507657747722648,0.473661427034500,0.441941738241593,0.412346222116530,0.384732629170287,0.358968235936574,0.334929207042592,0.312500000000000,0.291572809855253,0.272047051030039,0.253828873861324,0.236830713517250,0.220970869120796,0.206173111058265,0.192366314585143,0.179484117968287,0.167464603521296,0.156250000000000,0.145786404927626]
frequencies = np.array(frequencies)/fs # normalize
# covert to suitable frequency range
scalesRange = pywt.frequency2scale(wavelet, frequencies)
MaxCoefThres = 0.1  # Threshold for CWT coefficients (max)
MaxCoefThres2 = 0.05  # Threshold for CWT coefficients (mean)
MinWalkingPeriod = 5  # Minimum walking period in seconds

# Defin parameters for feature extraction
window = 10  # Window size in seconds
lyingOrient = 85  # Lying orientation in degrees
staMethod = 5  # Statistics method
WalkingPeriod = np.arange(8, 19)  # CWT periods for walking activity

#-----------------------------------------------------------------------------------
# Helper functions
# Define staset function
def staset(V, ax, method):
    if method == 1:
        mean_val = np.mean(V, axis=ax)
        std_val = np.std(V, axis=ax)
        cv_val = (std_val / mean_val) * 100 if mean_val != 0 else np.nan
        median_val = np.median(V, axis=ax)
        skew_val = abs(skew(V, axis=ax))
        kurt_val = kurtosis(V, axis=ax)
        staout = [mean_val, std_val, cv_val, median_val, skew_val, kurt_val]

    elif method == 5:
        mean_val = np.mean(V, axis=ax)
        std_val = np.std(V, axis=ax)
        cv_val = (std_val / mean_val) * 100 if mean_val != 0 else np.nan
        median_val = np.median(V, axis=ax)
        perc_25 = np.percentile(V, 25, axis=ax)
        perc_75 = np.percentile(V, 75, axis=ax)
        staout = [mean_val, std_val, cv_val, median_val, perc_25, perc_75]

    elif method==6:
        mean_val = np.mean(V, axis=ax)
        std_val = np.std(V, axis=ax)
        cv_val = (std_val / mean_val) * 100 if mean_val != 0 else np.nan
        median_val = np.median(V, axis=ax)
        staout = [mean_val, std_val, cv_val, median_val]
    return np.array(staout)

def find_stationary_periods(walking_coef, max_coef_thres, max_coef_thres2):
    stationary = (np.max(np.abs(walking_coef), axis=0) > max_coef_thres) & \
                 (np.mean(np.abs(walking_coef), axis=0) > max_coef_thres2)
    return stationary

#------------------------------------Step1: data preprocessing----------------------------------------

# Preallocate data structures
filtMagdata = {}
filtxyzdata = {}
cwtCoef = {}
Labels = {}

# Loop through participants to read and preprocess data
for i in range(len(ParID)):
    # Select filename based on participant ID
    if i < 18:
        filename = f"{ParID[i]}.csv"
        pathnamesel = pathname70
    elif ParID[i] < 10:
        filename = f"S00{ParID[i]}.csv"
        pathnamesel = pathname
    else:
        filename = f"S0{ParID[i]}.csv"
        pathnamesel = pathname
    
    # Read data
    data = pd.read_csv(f"{pathnamesel}/{filename}")
    
    # Filter magnitude and XYZ data
    mag_data = np.linalg.norm(data.iloc[:, 1:4], axis=1)
    filtdata = filtfilt(LPa, LPb, mag_data)
    filtdataxyz = filtfilt(LPa, LPb, data.iloc[:, 1:4].values, axis=0)

    # Perform CWT (Continuous Wavelet Transform) using PyWavelets
    print(f'Processing CWT for subject {i + 1}')
    coef, freqs = pywt.cwt(filtdata - 1, scales=scalesRange, wavelet=wavelet, sampling_period=1/fs)
    #scales=np.arange(1, 128),
    cwtCoef[f'Subject{i + 1}'] = coef
    filtMagdata[f'Subject{i + 1}'] = filtdata
    filtxyzdata[f'Subject{i + 1}'] = filtdataxyz
    Labels[f'Subject{i + 1}'] = data.iloc[:, -1].values
    print(f'Finished subject {i + 1}')

    # Save all objects into a single file
    processed_data = {
    'cwtCoef': cwtCoef,
    'filtMagdata': filtMagdata,
    'filtxyzdata': filtxyzdata,
    'Labels': Labels
    }
   #joblib.dump(processed_data, os.path.join('features', 'processed_data.pkl'))
   #print("Data saved as processed_data.pkl")
   

# -----------------------------------STEP2: feature extraction---------------------------------------------

RealLabel = []
ParticipantLabel = []
f = []
# Loop through each participant to extract features
for i in range(len(ParID)):       
    wceof = cwtCoef[f'Subject{i + 1}']
    datause = filtMagdata[f'Subject{i + 1}']
    datausexyz = filtxyzdata[f'Subject{i + 1}']
    
    # Estimate body orientation
    BodyOrientEst = np.abs(np.degrees(np.arctan(np.sqrt(np.sum(datausexyz[:, [1, 2]] ** 2, axis=1)) / np.abs(datausexyz[:, 0]))))
    ChangeOrient = np.abs(np.diff(np.insert(BodyOrientEst, 0, 0)))
    
    # Find stationary periods using CWT coefficients
    stationary = find_stationary_periods(wceof[WalkingPeriod, :], MaxCoefThres, MaxCoefThres2)
    
    # Identify edges of stationary periods and remove short sequences
    d_data = np.diff(np.concatenate([[0], stationary, [0]]))
    startIdx = np.where(d_data == 1)[0]
    endIdx = np.where(d_data == -1)[0] - 1
    lengths = endIdx - startIdx + 1
    
    # Remove short walking periods
    shortSeqIdx = np.where(lengths < MinWalkingPeriod * fs)
    for j in shortSeqIdx[0]:
        stationary[startIdx[j]:endIdx[j]] = 0

    # Window-based feature extraction
    idx = np.arange(0, len(datause), window * fs)
    if idx[-1] != len(datause):
        idx = np.append(idx, len(datause))

    # Extract features from each window
    for j in range(len(idx) - 1):
        timerange = np.arange(idx[j], idx[j + 1])
        if j == 0:
            ChangeOrient1 = ChangeOrient[timerange[1:]]
        else:
            ChangeOrient1 = ChangeOrient[timerange]

        BodyOrientEstuse = BodyOrientEst[timerange]

        absCoef = np.abs(wceof[:, timerange])
        absCoefWalking = np.abs(wceof[WalkingPeriod, :][:, timerange])

        # Feature extraction using peaks in the walking coefficient
        peaks1Loc, _ = find_peaks(np.mean(absCoefWalking, axis=0), height=0.05, prominence=0.01)
        peaks1=np.mean(absCoefWalking, axis=0)[peaks1Loc]
        if len(peaks1) == 0:
            staout1 = [0] * 9
        else:
            staout1 = []
            # Extend with elements from staset(peaks1, 0, staMethod)
            staout1.extend(staset(peaks1, 0, 6).flatten())
            # Extend with elements from staset(peaks1Loc, 0, staMethod)
            staout1.extend(staset(peaks1Loc, 0, 6).flatten())
            # Append the scalar value len(peaks1)
            staout1.append(len(peaks1))
            
        maxv_pr = np.max(np.max(absCoef, axis=1))
        maxv_wr = np.max(np.max(absCoefWalking, axis=1))
        maxvl_pr = np.argmax(np.max(absCoef, axis=1))
        maxvl_wr = np.argmax(np.max(absCoefWalking, axis=1))
        
        # Aggregate features using staset
        staout = np.concatenate([
            staset(datause[timerange], 0, staMethod),
            staset(datausexyz[timerange, 0], 0, staMethod),
            staset(np.max(absCoef, axis=0), 0, staMethod),
            staset(np.mean(absCoef, axis=0), 0, staMethod),
            staset(np.max(absCoefWalking, axis=0), 0, staMethod),
            staset(np.mean(absCoefWalking, axis=0), 0, staMethod),
            staset(np.mean(absCoefWalking, axis=0)/np.mean(absCoef, axis=0), 0, staMethod),
            staout1,
            np.array([1/frequencies[maxvl_pr], maxv_pr, 1/frequencies[maxvl_wr], maxv_wr]),
            staset(ChangeOrient1, 0, staMethod),
            staset(BodyOrientEstuse, 0, staMethod),
            np.array([np.max(BodyOrientEstuse) - np.min(BodyOrientEstuse), 
                      np.mean(BodyOrientEstuse) - lyingOrient, np.sum(stationary[timerange])])
        ])

        # Store labels and features
        labelint = mode(Labels[f'Subject{i + 1}'][timerange])
        RealLabel.append(labelint[0])
        ParticipantLabel.append(i + 1)
        f.append(staout)
    print(f'Finished extracting features for subject {i + 1}')
    

#----------------------------------- Step3: save features used for train model-------------------------------------
    
# Ensure features and labels are in the correct format
f_array = np.array(f)  # Convert list of features to a NumPy array
real_labels = np.array(RealLabel).flatten()  # Flatten real labels

# Initialize WalkingSSLO as a string array with the same shape as real_labels
WalkingSSLO = np.empty_like(real_labels, dtype=object)
# Apply conditions and assign labels
WalkingSSLO[np.isin(real_labels, [1, 4, 5])] = "Walking"
WalkingSSLO[np.isin(real_labels, [3, 6])] = "Standing"
WalkingSSLO[real_labels == 7] = "Sitting"
WalkingSSLO[real_labels == 8] = "Lying"
WalkingSSLO[np.isin(real_labels, [2, 13, 14, 130, 140])] = "Other"

participant_labels = np.array(ParticipantLabel).flatten()  # Flatten participant labels
# Create a DataFrame for the features
feature_names = [
    "AccMag_mean", "AccMag_std", "AccMag_var", "AccMag_median",
    "AccMag_25perc", "AccMag_75perc", "StandAxis_mean", "StandAxis_std",
    "StandAxis_var", "StandAxis_median", "StandAxis_25perc",
    "StandAxis_75perc", "ScMaxPAllf_mean", "ScMaxPAllf_std", "ScMaxPAllf_var",
    "ScMaxPAllf_median", "ScMaxPAllf_25perc", "ScMaxPAllf_75perc",
    "ScMeanPAllf_mean", "ScMeanPAllf_std", "ScMeanPAllf_var", "ScMeanPAllf_median",
    "ScMeanPAllf_25perc", "ScMeanPAllf_75perc", "ScMaxPWalkf_mean",
    "ScMaxPWalkf_std", "ScMaxPWalkf_var", "ScMaxPWalkf_median", "ScMaxPWalkf_25perc",
    "ScMaxPWalkf_75perc", "ScMeanPWalkf_mean", "ScMeanPWalkf_std", "ScMeanPWalkf_var",
    "ScMeanPWalkf_median", "ScMeanPWalkf_25perc", "ScMeanPWalkf_75perc",
    "ScRatioP_mean", "ScRatioP_std", "ScRatioP_var", "ScRatioP_median", "ScRatioP_25perc",
    "ScRatioP_75perc", "CoefPeaks_mean", "CoefPeaks_std", "CoefPeaks_var",
    "CoefPeaks_median", "CoefPeakWidths_mean", "CoefPeakWidths_std", "CoefPeakWidths_var",
    "CoefPeakWidths_median", "CoefPeaks_", "DomfAll_", "TiAveDomfAll_", "DomfWalk_", "TiAveDomfWalk_",
    "OrientChange_mean", "OrientChange_std", "OrientChange_var", "OrientChange_median",
    "OrientChange_25perc", "OrientChange_75perc", "Orient_mean", "Orient_std", "Orient_var",
    "Orient_median", "Orient_25perc", "Orient_75perc", "MaxOrientChange_", "MeanOrientDiffLying_", "Perctime_"
]
features_df = pd.DataFrame(f_array, columns=feature_names)

# Add labels and participant information
features_df['Labels'] = WalkingSSLO  # Primary labels
# Update WalkingSSLO
WalkingSSLO[np.isin(WalkingSSLO, ['Standing', 'Sitting'])] = 'Standing/Sitting'
features_df['Labels2'] = WalkingSSLO  # Primary labels
features_df['ParLabel'] = participant_labels  # Participant-specific labels

# Save the DataFrame to a CSV file
features_df.to_csv(os.path.join('features', 'processed_features_python.csv'), index=False)
print("Features saved to 'features\processed_features_python.csv'")



    