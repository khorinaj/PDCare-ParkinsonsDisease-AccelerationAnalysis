"""
Script Name: feature_selection_evaluation.py

Description:
This script evaluates the performance of a Random Forest classifier using features selected based on their importance scores.
It performs Leave-One-Participant-Out Cross-Validation (LOPOCV) and computes various metrics (e.g., precision, recall, F1-score) 
per activity, macro-averaged, and weighted-averaged. It also visualizes the confusion matrix and saves feature importance scores.

Inputs:
- Pre-trained models (`my_model4class.pkl` and `my_model4classfs36_pt2.pkl`). Change to your own model when required
- Dataset: `features/processed_features_matlab.csv` containing features, labels, and participant information.

Outputs:
- Confusion matrix as a percentage heatmap and a standard confusion matrix.
- Metrics summarized by activity and averaged across folds.
- Feature importance scores saved as an Excel file.

Author: Yanke Sun
Date: 31/10/2024
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import LeaveOneGroupOut,cross_val_predict
from sklearn.metrics import classification_report, f1_score,confusion_matrix,ConfusionMatrixDisplay,accuracy_score
import joblib
from joblib import Parallel, delayed
from matplotlib import pyplot as plt
import seaborn as sns  # Make sure to import seaborn

results_folder = 'results'
model1='my_model4class.pkl' # full features 
model2='my_model4classfs36_pt2.pkl'# selected features

# Load the model for feature selection
loaded_model = joblib.load(os.path.join(results_folder,model1))
best_model = loaded_model.best_estimator_

# Load data
dataset = pd.read_csv(os.path.join('features', 'processed_features_matlab.csv'))
fMat=dataset.iloc[:, :70]
classLabels=dataset[['Labels2']].values.flatten()
ParLabel = dataset[['ParLabel']].values.flatten() # Convert Series to a 1-D NumPy array

# LOPOCV
logo = LeaveOneGroupOut()
logo.get_n_splits(fMat, classLabels, ParLabel)
logo.get_n_splits(groups=ParLabel)  # 'groups' is always required


# Extract feature importances from the trained Random Forest model
feature_importances = best_model.feature_importances_

# Create a mask for selecting features based on the threshold
#important_features_mask = feature_importances > feature_importances.mean()
important_features_mask = feature_importances > 0.005
print('Number of features selected: ',sum(important_features_mask))

# save feature importance score 
file_path ='featureImportance.xlsx'
# Convert y_pred to a pandas DataFrame
y_pred_df = pd.DataFrame(feature_importances, columns=['feature_importances'])
y_pred_df.to_excel(os.path.join('results', file_path), index=False)

# Select important features
f =fMat.loc[:, important_features_mask]
# Load the model grid serach on feature selected

loaded_model = joblib.load(os.path.join(results_folder, model2))
best_model = loaded_model.best_estimator_
# Output the best parameters and estimator
print("Best parameters found: ", loaded_model.best_params_)
print("Best estimator: ", best_model)

# Perform LOOCV and collect predictions using the best estimator
y_pred = cross_val_predict(best_model, f, classLabels, cv=logo, groups=ParLabel, n_jobs=-1)
# Compute the confusion matrix
conf_matrix = confusion_matrix(classLabels, y_pred)

# Plot the confusion matrix
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix, display_labels=best_model.classes_)
disp.plot(cmap=plt.cm.Blues)
plt.savefig(os.path.join('results', 'CM4classfs36_pt2'), dpi=300, bbox_inches='tight')

# Convert to percentages
cm_percentage = conf_matrix.astype('float') / conf_matrix.sum(axis=1)[:, np.newaxis] *100
# Create a DataFrame for better visualization with seaborn
cm_df = pd.DataFrame(cm_percentage, index=best_model.classes_, columns=best_model.classes_)
def custom_annotation(data):
    annotations = data.copy()
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            annotations.iloc[i, j] = f"{data.iloc[i, j]:.1f}%"
    return annotations

annotations = custom_annotation(cm_df)

# Plotting and save the confusion matrix
plt.figure(figsize=(10, 7))
sns.set(font_scale=1.2)
sns.heatmap(cm_df, annot=annotations.values, fmt='', cmap='Blues', annot_kws={"size": 16})
plt.ylabel('True label')
plt.xlabel('Predicted label')
plt.savefig(os.path.join('results', 'CM4classfs36_pt2_percentage.png'), dpi=300, bbox_inches='tight')


# helper functions
def train_and_evaluate(train_index, test_index):
    X_train, y_train = f.iloc[train_index], classLabels[train_index]
    X_test, y_test = f.iloc[test_index], classLabels[test_index]
    
    best_model.fit(X_train, y_train)
    y_pred = best_model.predict(X_test)
    
    report = classification_report(y_test, y_pred, output_dict=True, zero_division=0)
    
    fold_metrics = {
        'precision': {},
        'recall': {},
        'f1-score': {},
        'accuracy': report['accuracy'],
        'macro_avg': report['macro avg'],
        'weighted_avg': report['weighted avg']
    }
    
    for key in ['Lying', 'Sit\\Stand', 'Walking', 'Other']:
        if key in report:
            fold_metrics['precision'][key] = report[key]['precision']
            fold_metrics['recall'][key] = report[key]['recall']
            fold_metrics['f1-score'][key] = report[key]['f1-score']
        else:
            fold_metrics['precision'][key] = np.nan
            fold_metrics['recall'][key] = np.nan
            fold_metrics['f1-score'][key] = np.nan
    
    participant_id = ParLabel[test_index][0]
    
    return fold_metrics, participant_id

results = Parallel(n_jobs=-1)(delayed(train_and_evaluate)(train_index, test_index) 
                              for train_index, test_index in logo.split(f, groups=ParLabel))

# Initialize dictionaries to store metrics
precision = {key: [] for key in ['Lying', 'Sit\\Stand', 'Walking', 'Other']}
recall = {key: [] for key in ['Lying', 'Sit\\Stand', 'Walking', 'Other']}
f1 = {key: [] for key in ['Lying', 'Sit\\Stand', 'Walking', 'Other']}
macro_avg = {'precision': [], 'recall': [], 'f1-score': []}
weighted_avg = {'precision': [], 'recall': [], 'f1-score': []}
acc = []

participant_metrics = {key: {} for key in ['Lying', 'Sit\\Stand', 'Walking', 'Other']}

for fold_metrics, participant_id in results:
    for key in ['Lying', 'Sit\\Stand', 'Walking', 'Other']:
        precision[key].append(fold_metrics['precision'][key])
        recall[key].append(fold_metrics['recall'][key])
        f1[key].append(fold_metrics['f1-score'][key])
        
        if participant_id not in participant_metrics[key]:
            participant_metrics[key][participant_id] = {'precision': [], 'recall': [], 'f1-score': []}
        participant_metrics[key][participant_id]['precision'].append(fold_metrics['precision'][key])
        participant_metrics[key][participant_id]['recall'].append(fold_metrics['recall'][key])
        participant_metrics[key][participant_id]['f1-score'].append(fold_metrics['f1-score'][key])
    
    acc.append(fold_metrics['accuracy'])
    for key in ['precision', 'recall', 'f1-score']:
        macro_avg[key].append(fold_metrics['macro_avg'][key])
        weighted_avg[key].append(fold_metrics['weighted_avg'][key])

# Calculate MetricPerAct
listperAct1 = {key: [] for key in ['precision', 'recall', 'f1-score']}
listperAct2 = {key: [] for key in ['precision', 'recall', 'f1-score']}
for key in precision.keys():
    listperAct1['precision'].append(np.nanmean(precision[key]))
    listperAct1['recall'].append(np.nanmean(recall[key]))
    listperAct1['f1-score'].append(np.nanmean(f1[key]))
    listperAct2['precision'].append(np.nanstd(precision[key]))
    listperAct2['recall'].append(np.nanstd(recall[key]))
    listperAct2['f1-score'].append(np.nanstd(f1[key]))

MetricPerAct = pd.concat([pd.DataFrame(listperAct1), pd.DataFrame(listperAct2)], axis=1)
MetricPerAct.index = ['Lying', 'Sit\\Stand', 'Walking', 'Other']
MetricPerAct.columns = ['precision_mean', 'recall_mean', 'f1-score_mean', 
                        'precision_std', 'recall_std', 'f1-score_std']

# Calculate MetricPerAct2
averaged_list1 = [np.mean(acc)] + [np.mean(macro_avg[key]) for key in macro_avg]
averaged_list2 = [np.std(acc)] + [np.std(macro_avg[key]) for key in macro_avg]
MetricPerAct2 = pd.concat([pd.DataFrame(averaged_list1), pd.DataFrame(averaged_list2)], axis=1)
MetricPerAct2.index = ['accuracy', 'precision', 'recall', 'f1-score']
MetricPerAct2.columns = ['mean', 'std']

# Calculate MetricPerAct3
weighted_list1 = [np.mean(weighted_avg[key]) for key in weighted_avg]
weighted_list2 = [np.std(weighted_avg[key]) for key in weighted_avg]
MetricPerAct3 = pd.concat([pd.DataFrame(weighted_list1), pd.DataFrame(weighted_list2)], axis=1)
MetricPerAct3.index = ['precision', 'recall', 'f1-score']
MetricPerAct3.columns = ['mean', 'std']

# Print results
print("MetricPerAct (Fold-averaged metrics per activity):")
print(MetricPerAct)
print("\nMetricPerAct2 (Overall accuracy and macro-averaged metrics):")
print(MetricPerAct2)
print("\nMetricPerAct3 (Weighted-averaged metrics):")
print(MetricPerAct3)