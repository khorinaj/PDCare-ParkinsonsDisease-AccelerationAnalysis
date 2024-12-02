"""
Script Name: train_model.py

Description:
This script trains and evaluates a Random Forest Classifier using Leave-One-Participant-Out Cross-Validation (LOPOCV).
The features used in this project are extracted using a MATLAB script and saved in 
'processed_features_matlab.csv'. While an equivalent Python script for feature extraction exists, 
there may be slight differences in the extracted features due to variations in continuous wavelet 
transform (CWT) functions and wavelets used in MATLAB and Python.

Key Steps:
1. Load pre-extracted features, class labels, and participant IDs from a CSV file.
2. Perform LOPOCV to train and validate the Random Forest model.
3. Optimize hyperparameters using Randomized Search in the initial round.
4. Conduct feature selection based on the importance scores from the first round.
5. Extract selected features, fine-tune the Random Forest model using Grid Search, and evaluate the final model.

Inputs:
- Feature datasets:
  - `features/processed_features_matlab.csv`: Features extracted using the MATLAB script.
  OR - `features/processed_features_python.csv`: Features extracted using the Python script 
     (may have minor differences due to implementation nuances).

Outputs:
- Serialized model files:
  - my_modelNew.pkl` (first model with full features).
  - my_modelNew_fs.pkl` (final fine-tuned model with selected features).
- Cross-validation predictions for the evaluation phase for two models.
- Confusion matrix visualization saved as an image for two models.

Additional Notes:
This is the final optimized model trained and saved by this study and can be used directly.
- **First Model (Full Feature Set):** The original model trained with the full feature set is saved 
  as `my_model4class.pkl`.
- **Second Model (Selected Features):** The model trained with the selected feature set and fine-tuned 
  parameters is saved as `my_model4classfs36.pkl`.
Author: Yanke Sun
Date: 31/10/2024

"""

import numpy as np
import pandas as pd
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_curve, auc,accuracy_score,confusion_matrix,ConfusionMatrixDisplay
from sklearn.model_selection import LeaveOneGroupOut,GridSearchCV,RandomizedSearchCV,cross_val_predict
from matplotlib import pyplot as plt
import joblib

# Load data
dataset = pd.read_csv(os.path.join('features', 'features_labels_participants.csv'))
fMat=dataset.iloc[:, :70]
classLabels=dataset[['Labels2']].values.flatten()
ParLabel = dataset[['ParLabel']].values.flatten() # Convert Series to a 1-D NumPy array

# LOPOCV 
logo = LeaveOneGroupOut()
logo.get_n_splits(fMat, classLabels, ParLabel)
logo.get_n_splits(groups=ParLabel)  # 'groups' is always required

model = RandomForestClassifier(random_state=42,n_jobs=-1)

# Define the parameter grid
param_grid = {
    'n_estimators': [250,350,450,550,650],
    'max_features': ['sqrt', 'log2'],
    'max_depth': [6,12,18,24],
    'min_samples_split': [2,4,6,8],
    'min_samples_leaf': [1, 2, 4],
    'bootstrap': [True, False],
    'criterion': ['gini', 'entropy']
}

f =fMat.loc[:, :]
# Create RandomizedSearchCV with the defined param_dist
grid_search = RandomizedSearchCV(estimator=model, param_distributions=param_grid, n_iter=200, cv=logo, scoring='f1_macro', n_jobs=-1, verbose=10,random_state=42)
#grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=logo, scoring='f1_macro', n_jobs=-1, verbose=10)
# Fit the model
grid_search.fit(f,classLabels,groups=ParLabel)

results_folder = 'results'
# Ensure the folder exists
os.makedirs(results_folder, exist_ok=True)
joblib.dump(grid_search, os.path.join(results_folder, 'my_modelNew.pkl'))

# Output the best parameters and estimator
print("Best parameters found: ", grid_search.best_params_)
print("Best estimator: ", grid_search.best_estimator_)

# Extract the best estimator from the grid search
best_model = grid_search.best_estimator_

# Perform LOOCV and collect predictions using the best estimator
y_pred = cross_val_predict(best_model, fMat, classLabels, cv=logo, groups=ParLabel, n_jobs=-1)

# Compute the confusion matrix
conf_matrix = confusion_matrix(classLabels, y_pred)
# Plot the confusion matrix
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix, display_labels=best_model.classes_)
disp.plot(cmap=plt.cm.Blues)
plt.savefig(os.path.join(results_folder, 'CM_my_modelNew'), dpi=300, bbox_inches='tight')


#------------------------------ second round after feature selection --------------------------------------------
# Extract feature importances from the trained Random Forest model
feature_importances = best_model.feature_importances_

# Create a mask for selecting features based on the threshold
important_features_mask = feature_importances > 0.005
print('Number of features selected: ',sum(important_features_mask))

# Select important features
f =fMat.loc[:, important_features_mask]

# Load the model grid serach on feature selected
param_grid = {
    'min_samples_split': [2,3,4,5],
    'max_depth': [8 , 12,  16,  20, 24, 28, 32, 36, 40],
    'n_estimators': [450,500,550]
}

# Create RandomizedSearchCV with the defined param_dist
grid_search = GridSearchCV(estimator=model, param_grid=param_grid, cv=logo, scoring='f1_macro', n_jobs=-1, verbose=10)
# Fit the model
grid_search.fit(f,classLabels,groups=ParLabel)

results_folder = 'results'
# Ensure the folder exists
os.makedirs(results_folder, exist_ok=True)
joblib.dump(grid_search, os.path.join(results_folder, 'my_modelNew_fs.pkl'))

# Output the best parameters and estimator
print("Best parameters found: ", grid_search.best_params_)
print("Best estimator: ", grid_search.best_estimator_)

# Extract the best estimator from the grid search
best_model = grid_search.best_estimator_

# Perform LOOCV and collect predictions using the best estimator
y_pred = cross_val_predict(best_model, fMat, classLabels, cv=logo, groups=ParLabel, n_jobs=-1)

# Compute the confusion matrix
conf_matrix = confusion_matrix(classLabels, y_pred)
# Plot the confusion matrix
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix, display_labels=best_model.classes_)
disp.plot(cmap=plt.cm.Blues)
plt.savefig(os.path.join(results_folder, 'CM_my_modelNew_fs'), dpi=300, bbox_inches='tight')