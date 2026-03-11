######################################################################
## This repository holds the codes for the research paper:          ##
##                                                                  ##
## Paper Title:                                                     ##
##   AI-Driven Discovery and Single-Cell Validation of a            ##
##   Conserved Epithelial Signature in HPV-Positive Malignancies    ##
##                                                                  ##
## File Author:                                                     ##
##   Naisarg Patel                                                  ##
##                                                                  ##
##                                                                  ##
## License: GNU General Public License v3.0 (GPL-3.0)               ##
## Contact: naisargbpatel14<at>gmail<dot>com                        ##
######################################################################   


import pandas as pd
import numpy as np
import lightgbm as lgb
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm
import os


def train_lgbm_single_feature(X_feature, y, feature_name, dataset_name):
    """
    Train LightGBM on a single feature and return metrics
    """
    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(
        X_feature.values.reshape(-1, 1), y, test_size=0.3, random_state=42, stratify=y
    )
    
    # Initialize LightGBM
    lgbm = lgb.LGBMClassifier(
        n_estimators=100,
        random_state=42,
        verbose=-1,  # Suppress LightGBM warnings
        force_col_wise=True,
        objective='binary' if len(np.unique(y)) == 2 else 'multiclass'
    )
    
    # Train the model
    lgbm.fit(X_train, y_train)
    
    # Make predictions
    y_pred = lgbm.predict(X_test)
    y_pred_proba = lgbm.predict_proba(X_test)[:, 1] if len(np.unique(y)) == 2 else None
    
    # Calculate metrics
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred, average='weighted')
    recall = recall_score(y_test, y_pred, average='weighted')
    f1 = f1_score(y_test, y_pred, average='weighted')
    
    # ROC AUC for binary classification only
    roc_auc = None
    if len(np.unique(y)) == 2 and y_pred_proba is not None:
        try:
            roc_auc = roc_auc_score(y_test, y_pred_proba)
        except Exception as e:
            print(f"  Warning: Could not calculate ROC AUC for feature {feature_name} in dataset {dataset_name}. Error: {str(e)}")
            roc_auc = None
    
    # Cross-validation score
    cv_scores = cross_val_score(lgbm, X_feature.values.reshape(-1, 1), y, 
                               cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42), 
                               scoring='accuracy')
    cv_mean = cv_scores.mean()
    cv_std = cv_scores.std()
    
    return {
        'Dataset': dataset_name,
        'Feature': feature_name,
        'Accuracy': accuracy,
        'Precision': precision,
        'Recall': recall,
        'F1_Score': f1,
        'ROC_AUC': roc_auc,
        'CV_Mean_Accuracy': cv_mean,
        'CV_Std_Accuracy': cv_std,
        'Train_Size': len(X_train),
        'Test_Size': len(X_test),
        'Feature_Importance': lgbm.feature_importances_[0] if len(lgbm.feature_importances_) > 0 else None
    }

def process_dataset(file_path, dataset_name):
    """
    Process a dataset and train LGBM on each feature
    """
    print(f"\nProcessing {dataset_name} dataset...")
    
    # Read the CSV file
    df = pd.read_csv(file_path)
    print(f"Dataset shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    
    # Prepare target variable
    le = LabelEncoder()
    y = le.fit_transform(df['tissue_type'])
    print(f"Target classes: {le.classes_}")
    print(f"Target distribution: {np.bincount(y)}")
    
    # Get feature columns (exclude SampleID and tissue_type)
    feature_columns = [col for col in df.columns if col not in ['SampleID', 'tissue_type']]
    print(f"Number of features: {len(feature_columns)}")
    
    results = []
    
    # Train LGBM on each feature individually
    for feature in tqdm(feature_columns, desc=f"Processing {dataset_name} features", unit="feature"):
        # Get the feature data
        X_feature = df[feature]
        
        # Check for missing values
        if X_feature.isnull().sum() > 0:
            tqdm.write(f"  Warning: Feature {feature} has {X_feature.isnull().sum()} missing values. Filling with median.")
            X_feature = X_feature.fillna(X_feature.median())
        
        # Train LGBM and get metrics
        try:
            metrics = train_lgbm_single_feature(X_feature, y, feature, dataset_name)
            results.append(metrics)
        except Exception as e:
            tqdm.write(f"  Error training on feature {feature}: {str(e)}")
            continue
    
    return results

def main():
    """
    Main function to process both datasets
    """
    print("Starting LightGBM Single Feature Analysis...")
    
    # File paths
    datasets = {
        'CC': 'CC.csv',
        'HNC': 'HNC.csv'
    }
    
    all_results = []
    
    # Process each dataset with overall progress
    with tqdm(total=len(datasets), desc="Processing datasets", unit="dataset") as pbar:
        for dataset_name, file_path in datasets.items():
            try:
                results = process_dataset(file_path, dataset_name)
                all_results.extend(results)
                pbar.set_postfix({"Current": dataset_name, "Features processed": len(results)})
            except Exception as e:
                tqdm.write(f"Error processing {dataset_name}: {str(e)}")
                continue
            finally:
                pbar.update(1)
    
    # Convert results to DataFrame
    results_df = pd.DataFrame(all_results)
    
    if not results_df.empty:
        # Sort by accuracy (descending)
        results_df = results_df.sort_values('Accuracy', ascending=False)
        
        # Create comparison DataFrame with requested columns
        features = set(results_df['Feature'].unique())
        comparison_data = []
        
        for feature in features:
            feature_data = results_df[results_df['Feature'] == feature]
            
            # Get accuracy for each dataset
            cc_data = feature_data[feature_data['Dataset'] == 'CC']
            hnc_data = feature_data[feature_data['Dataset'] == 'HNC']
            
            cc_accuracy = cc_data['Accuracy'].iloc[0] if len(cc_data) > 0 else None
            hnc_accuracy = hnc_data['Accuracy'].iloc[0] if len(hnc_data) > 0 else None
            
            # Calculate average accuracy (only if both datasets have data)
            if cc_accuracy is not None and hnc_accuracy is not None:
                avg_accuracy = (cc_accuracy + hnc_accuracy) / 2
            else:
                avg_accuracy = cc_accuracy if cc_accuracy is not None else hnc_accuracy
            
            comparison_data.append({
                'Feature_Name': feature,
                'Accuracy_in_CC': cc_accuracy,
                'Accuracy_in_HNC': hnc_accuracy,
                'Average_Accuracy': avg_accuracy
            })
        
        # Create comparison DataFrame
        comparison_df = pd.DataFrame(comparison_data)
        
        # Calculate ranks
        if 'Accuracy_in_CC' in comparison_df.columns:
            comparison_df['Rank_in_CC'] = comparison_df['Accuracy_in_CC'].rank(method='dense', ascending=False, na_option='bottom').astype('Int64')
        
        if 'Accuracy_in_HNC' in comparison_df.columns:
            comparison_df['Rank_in_HNC'] = comparison_df['Accuracy_in_HNC'].rank(method='dense', ascending=False, na_option='bottom').astype('Int64')
        
        comparison_df['Rank_based_on_Average_Accuracy'] = comparison_df['Average_Accuracy'].rank(method='dense', ascending=False, na_option='bottom').astype('Int64')
        
        # Sort by average accuracy rank
        comparison_df = comparison_df.sort_values('Rank_based_on_Average_Accuracy')
        
        # Create output folder
        output_folder = 'LGBM_Results'
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
            print(f"Created output folder: {output_folder}")
        
        # Save comparison results to CSV
        output_file_combined = os.path.join(output_folder, 'lgbm_single_feature_results_combined.csv')
        comparison_df.to_csv(output_file_combined, index=False)
        print(f"\nCombined comparison results saved to: {output_file_combined}")
        
        # Save separate files for each dataset (detailed results)
        output_files = {}
        for dataset in results_df['Dataset'].unique():
            dataset_results = results_df[results_df['Dataset'] == dataset].copy()
            dataset_results = dataset_results.sort_values('Accuracy', ascending=False)
            output_file = os.path.join(output_folder, f'lgbm_single_feature_results_{dataset}.csv')
            dataset_results.to_csv(output_file, index=False)
            output_files[dataset] = output_file
            print(f"{dataset} dataset results saved to: {output_file}")
        
        # Display summary statistics
        print("\n" + "="*80)
        print("SUMMARY STATISTICS")
        print("="*80)
        
        for dataset in results_df['Dataset'].unique():
            dataset_results = results_df[results_df['Dataset'] == dataset]
            print(f"\n{dataset} Dataset:")
            print(f"  Number of features: {len(dataset_results)}")
            print(f"  Best accuracy: {dataset_results['Accuracy'].max():.4f}")
            print(f"  Mean accuracy: {dataset_results['Accuracy'].mean():.4f}")
            print(f"  Std accuracy: {dataset_results['Accuracy'].std():.4f}")
            print(f"  Best feature: {dataset_results.iloc[0]['Feature']}")
        
        # Display top 10 features overall
        print("\n" + "="*80)
        print("TOP 10 FEATURES (by Average Accuracy)")
        print("="*80)
        top_10 = comparison_df.head(10)
        for idx, row in top_10.iterrows():
            cc_acc = f"{row['Accuracy_in_CC']:.4f}" if pd.notna(row['Accuracy_in_CC']) else 'N/A'
            hnc_acc = f"{row['Accuracy_in_HNC']:.4f}" if pd.notna(row['Accuracy_in_HNC']) else 'N/A'
            avg_acc = f"{row['Average_Accuracy']:.4f}" if pd.notna(row['Average_Accuracy']) else 'N/A'
            print(f"Rank {int(row['Rank_based_on_Average_Accuracy']):2d} | {row['Feature_Name']:<15} | CC: {cc_acc} | HNC: {hnc_acc} | Avg: {avg_acc}")
        
        print(f"\nTotal features analyzed: {len(comparison_df)}")
        print(f"Combined comparison results saved to: {output_file_combined}")
        for dataset, file_path in output_files.items():
            print(f"{dataset} detailed results saved to: {file_path}")
        print(f"All results saved in folder: {output_folder}")
    else:
        print("No results to save!")

if __name__ == "__main__":
    main()
