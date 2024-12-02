# Parkinson’s Disease activity features Analysis

## Description
This project aims to leverage machine learning models to predict human activities using accelerometers worn on the trunk or wrist, analyze long-term activity patterns in Parkinson’s Disease (PD) patients, and explore correlations between activity features and clinical scales. The ultimate objective is to identify potential features for long-term monitoring of PD symptoms and progression.

This project combines:
1. Human Activity Recognition (HAR):
- Use machine learning models trained on publicly available datasets with trunk sensor to classify activities like walking, sitting/standing, lying, and other behaviors.
2. Parkinson’s Disease (PD) Analysis:
- Preprocess real-world, long-term trunk and wrist sensor data collected from PD patients and find correlation between activity related features (based on predictions of activity on ML model) and clinical scales

## Folder Structure
```plaintext
├── data/                      # Raw PD data collected in our project (publicly available, not included in this repository)
├── HARmodel/                 # HAR analysis using Python
│   ├── src/                   # Python scripts for HAR pipeline
│   ├── results/               # Results from HAR pipeline
│   ├── features/              # Extracted features and labels
│   ├── main.py                # main script
│   ├── README.md              # Documentation for HAR pipeline
│   ├── requirements.txt           # Python dependencies
├── PDanalysis/                   # PD analysis 
│   ├──CorrelationPlot.ipynb   # a Juypter notebook for plotting correlation map used in the paper
│   ├── clinicalScales/        # Clinical scale data
│   ├── dataProcessed/         # Processed trunk sensor data
│   ├── utils/                 # MATLAB utility scripts
│   ├── results/               # Analysis results for trunk data
│   ├── main.m                 # main script
│   ├── README.md              # Documentation for trunk sensor analysis

```

## Subprojects

### HARmodel (for trunk sensor)
Description: Human activity recognition using Python and trunk sensor with publicly available datasets, HARTH and HAR70+, available in [NTNU AI-Lab HARTH repository](https://github.com/ntnu-ai-lab/harth-ml-experiments).
#### Key Features:
- Extracts features and trains a Random Forest classifier for activity prediction (ModelBody).
- Feature selection based on importance score and evaluates performance using cross-validation.
- Prediction on external PD data 

#### Outputs:
- Trained machine learning model.
- Confusion matrices and feature importance visualizations.
- Predicted labels for PD data
**Details of running the project refer to HARmodel/README.md for detailed instructions.**

### PDanalysis
Description: PD data analysis using trunk sensors and wirst sensors.
#### Key Features:
* Extract labels in long-term trunk sensor data predicted by the model trained in GitHub_PD (ModelBody).
   or Extract labels in in long-term wrist sensor data predicted by publicly available tools from [OxWearables](https://github.com/OxWearables/ssl-wearables).
* Extract activity-related features (e.g., percentage of sleep, mean walking acceleration) for 24-hour or 48-hour intervals.
* Performs correlation analysis between activity features and clinical scales to explore their potential for symptom monitoring and disease progression tracking..

#### Outputs:
- Activity features for each window.
- Correlation matrices and visualizations.
**Details of running the project refer to MATLAB/README.md for detailed instructions.**


## License
This project is licensed under the MIT License.

## Author
Yanke Sun
Email: khorinaj@outlook.com

## Citation
