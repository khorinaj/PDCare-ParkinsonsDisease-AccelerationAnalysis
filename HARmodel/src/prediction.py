"""
Script Name: prediction.py

Description:
This script loads pre-trained Random Forest models
Read MATLAB pre-processed data files for PD long-term dialy life data
Predicts activity labels for specified PD patient IDs. 
Predictions are saved as Excel and CSV files in the same directory as the MATLAB files.

Inputs:
- MATLAB `.mat` files containing processed features (`ID*_processed.mat`).
- Pre-trained Random Forest model (e.g `my_model4classfs36_pt2.pkl`).

Outputs:
- Activity predictions saved as Excel files (`ID*_Label4classfs_pt2.xlsx`).

Author: Yanke Sun
Date: 31/10/2024
"""


import numpy as np
import pandas as pd
import joblib
import scipy.io
import os

ids=["23","27","05","16","06","07","09"];


results_folder = 'results'
model1 = 'my_model4class.pkl'  # Full feature model
model2 = 'my_model4classfs36_pt2.pkl'  # Selected features model

# Ensure the results folder exists
os.makedirs(results_folder, exist_ok=True)

# Load the first model for feature selection
loaded_model = joblib.load(os.path.join(results_folder, model1))
best_model = loaded_model.best_estimator_

# Extract feature importances
feature_importances = best_model.feature_importances_
# Create a mask for selecting features with importance > 0.005
important_features_mask = feature_importances > 0.005

# Load the second model on feature selected
loaded_model = joblib.load(os.path.join(results_folder, model2))
best_model = loaded_model.best_estimator_

# Define the base folder dynamically based on the current file's location
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Path to the dataProcessed folder in the MATLAB directory
directory = os.path.join(BASE_DIR, 'PDanalysis', 'dataProcessed_Body')
# directory = r'GitHub Paper\PDanalysis\datapProcessed_Body'

for id in ids:
    file_name='ID'+ id +'_processed.mat'
    file_name = os.path.join(directory, file_name)
    MatlabData= scipy.io.loadmat(file_name)
    Var = MatlabData['fPD']
    df = pd.DataFrame(Var)
    fPD1=df.iloc[:,:70]
    # Select important features
    fselected1 = fPD1.loc[:, important_features_mask]
    y_pred = best_model.predict(fselected1)
    # Define the file path
    file_path ='ID'+ id + 'Label4classfs_pt2.xlsx'
    file_path=os.path.join(directory, file_path)
    # Convert y_pred to a pandas DataFrame
    y_pred_df = pd.DataFrame(y_pred, columns=['Prediction'])
    # Write the DataFrame to a CSV file
    y_pred_df.to_excel(file_path, index=False)
    print(y_pred)           
    