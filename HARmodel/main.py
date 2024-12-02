
# ----------------------------------- Script Description -----------------------------------
"""
Main Script: main.py

This script serves as the entry point for the pipeline. It orchestrates the entire process of
data preprocessing, model training, evaluation, and prediction by sequentially running individual
scripts located in the `src` directory.

Workflow:
1. Data Preprocessing: Prepares the dataset by extracting features from the public HAR dataset.
2. Model Training: Trains the Random Forest model using preprocessed features.
3. Model Evaluation: Evaluates the trained model using LOPOCV and generates performance metrics.
4. Prediction: Uses the trained model to predict activity labels on a PD dataset.

Dependencies:
- Ensure all scripts (`preprocess_feature_extraction.py`, `train_model.py`, `model_evaluation.py`, `prediction.py`)
  are located in the `src` folder and properly configured.
- Python 3.x

Author: Yanke Sun
Date: 31/10/2024
"""
# main.py
import os

# ---------------------------------- Set Working Directory ---------------------------------
# Ensure the script runs from the correct directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# ----------------------------------- Step 1: Preprocessing ---------------------------------
print("Running data preprocessing and feature extractions...")
exec(open("src/preprocess_feature_extraction.py").read())
print("Data preprocessing and feature extractions completed.")

# ----------------------------------- Step 2: Model Training --------------------------------
print("Running model training...")
exec(open("src/train_model.py").read())
print("Model training completed.")

# ----------------------------------- Step 3: Evaluation ------------------------------------
print("Running model evaluation...")
exec(open("src/model_evaluation.py").read())
print("Model evaluation completed.")

# ----------------------------------- Step 4: Prediction ------------------------------------
print("Running prediction...")
exec(open("src/prediction.py").read())
print("Prediction completed.")

# -------------------------------------- Completion -----------------------------------------
print("Pipeline completed.")


