import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os, sys, subprocess, random, pickle, warnings, getopt, statistics
from matplotlib.pyplot import *
from pandas import read_csv
from time import time
from numpy import asarray, sqrt, argmax
from statistics import median
from itertools import *
import sv_dataframe_transform
from sklearn import svm, metrics, datasets, preprocessing
from sklearn.preprocessing import LabelEncoder, label_binarize, StandardScaler, OrdinalEncoder, OneHotEncoder 
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.metrics import mean_squared_error, median_absolute_error, roc_curve, auc, accuracy_score, precision_score, recall_score, precision_recall_curve, confusion_matrix
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer, TransformedTargetRegressor
from sklearn.pipeline import make_pipeline
from supervised.automl import AutoML 

# Ignore warnings
warnings.filterwarnings("ignore")

############################################################################################################################################################# 
# Get current working directory and parse command line arguments
wd = sys.path[0]
opts, args = getopt.getopt(sys.argv[1:], "s:o:k:")
for op, value in opts:
    if op == "-s":
        cohort_name = str(value)
    if op == "-o":
        outdir = str(value)
    if op == "-k":
        kfolds = int(value)
sv_type_vector = ['trs', 'nontrs']
# subprocess.call(["sudo", "mkdir", "-p", outdir, '/cyto_sv_ml'])

############################################################################################################################################################# 
# Read input data
sv_data = pd.read_csv(outdir + '/' + str(cohort_name) + '.sv.all.combine_all', sep='\t', header=0, index_col=None, keep_default_na=False)
print(sv_data.columns)

# Process data based on sv_type
for sv_type in sv_type_vector:
    if sv_type == 'trs':
        # Filter out certain sv_types and save the resulting data
        sv_data_01 = sv_data[~(sv_data['sv_type'].isin(['DEL', 'DUP', 'INV', 'INS']))].copy()  
        sv_data_01.to_csv(outdir + '/' + str(cohort_name) + '.sv.all.combine_all_trs_o', sep='\t', header=True, index=None)        
        sv_data_1 = sv_dataframe_transform.trs_sv_data_transform(sv_data_01)[0]
        sv_data_1.to_csv(outdir + '/' + str(cohort_name) + '.sv.all.combine_all_trs', sep='\t', header=True, index=None)
        sv_data_10 = sv_dataframe_transform.trs_sv_data_transform(sv_data_01)[1]
        sv_data_10.to_csv(outdir + '/' + str(cohort_name) + '.sv.all.combine_all_trs_all', sep='\t', header=True, index=None)        
    else:
        # Filter in certain sv_types and save the resulting data
        sv_data_01 = sv_data[sv_data['sv_type'].isin(['DEL', 'DUP', 'INV', 'INS'])].copy()     
        sv_data_01.to_csv(outdir + '/' + str(cohort_name) + '.sv.all.combine_all_nontrs_o', sep='\t', header=True, index=None)       
        sv_data_1 = sv_dataframe_transform.nontrs_sv_data_transform(sv_data_01)[0]
        sv_data_1.to_csv(outdir + '/' + str(cohort_name) + '.sv.all.combine_all_nontrs', sep='\t', header=True, index=None)        
        sv_data_10 = sv_dataframe_transform.nontrs_sv_data_transform(sv_data_01)[1]
        sv_data_10.to_csv(outdir + '/' + str(cohort_name) + '.sv.all.combine_all_nontrs_all', sep='\t', header=True, index=None)            
    
    # Generate and save summary plot
    sv_summary_plot = sv_dataframe_transform.sv_data_summary_plot(sv_data_1)    
    sv_summary_plot.savefig(outdir + '/cyto_sv_ml/' + str(cohort_name) + '_' + sv_type + '_' + 'sv_summary_plot.pdf', transparent=True)  

    print("# Load SV data matrix with label")
    sv_data12_2 = sv_data_1[sv_data_1['label'] != 0].copy()
    sv_data12_2 = sv_data12_2.drop(['sv_type'], axis=1)
    sv_data12_0 = sv_data_1[sv_data_1['label'] == 0].copy()
    sv_data12_2.index = range(sv_data12_2.shape[0])
    y_false_sv = np.array(sv_data12_2[sv_data12_2['label'] == -1].index)
    y_germline_sv = np.array(sv_data12_2[sv_data12_2['label'] == 1].index)
    y_somatic_sv = np.array(sv_data12_2[sv_data12_2['label'] == 2].index)

    print("# Define categorical and continuous features")
    enc = preprocessing.LabelEncoder()
    categorical_varaibles = []  # Specify categorical variables here if any
    categorical_columns = [i for i in sv_data12_2.columns if i in categorical_varaibles]
    numerical_columns = [i for i in sv_data12_2.columns if i not in categorical_varaibles + ['label']]
    
    print("# Set up the pipeline for data transformation and modeling")
    preprocessor = make_column_transformer(
        (OneHotEncoder(drop="if_binary"), categorical_columns),
        (StandardScaler(), numerical_columns),
        remainder="passthrough"
    )
    model1 = make_pipeline(preprocessor)
                   
############################################################################################################################################################# 
                    
    print("# Tune sub-oversampling and split training/testing/validation") 
    y = sv_data12_2.label.values
    X = sv_data12_2.drop(["label"], axis=1)
    
    # Open the text file for saving model performance metrics
    sv_ml_metrics_ts_file = open(outdir + '/cyto_sv_ml/' + str(cohort_name) + '_' + sv_type + '_' + 'sv_ml_metrics_sub.csv', 'w')  
    sv_ml_metrics_ts = {'mic_pre': [], 'mac_pre': [], 'mic_rec': [], 'mac_rec': []}  
    
    label_avg = median([round(len(y_somatic_sv) * 0.9), round(len(y_germline_sv) * 0.9), round(len(y_false_sv) * 0.9)])
    
    for k in range(kfolds):
        np.random.seed(k * 2)
        
        # Subsample each class to balance the dataset
        if label_avg <= round(len(y_false_sv) * 0.9):
            false_sv_idx = np.random.choice(np.random.choice(y_false_sv, round(len(y_false_sv) * 0.9), replace=False), label_avg, replace=False)    
        else:
            false_sv_idx = np.random.choice(np.random.choice(y_false_sv, round(len(y_false_sv) * 0.9), replace=False), label_avg, replace=True)         
        
        if label_avg <= round(len(y_germline_sv) * 0.9):
            germline_sv_idx = np.random.choice(np.random.choice(y_germline_sv, round(len(y_germline_sv) * 0.9), replace=False), label_avg, replace=False)
        else:
            germline_sv_idx = np.random.choice(np.random.choice(y_germline_sv, round(len(y_germline_sv) * 0.9), replace=False), label_avg, replace=True)                   
        
        if label_avg <= round(len(y_somatic_sv) * 0.9):
            somatic_sv_idx = np.random.choice(np.random.choice(y_somatic_sv, round(len(y_somatic_sv) * 0.9), replace=False), label_avg, replace=False)
        else:                
            somatic_sv_idx = np.random.choice(np.random.choice(y_somatic_sv, round(len(y_somatic_sv) * 0.9), replace=False), label_avg, replace=True)
        
        idx = np.r_[germline_sv_idx, false_sv_idx, somatic_sv_idx]

        y12 = sv_data12_2.iloc[idx, :].label.values
        X12 = sv_data12_2.iloc[idx, :].drop(["label"], axis=1)
        y12_2 = sv_data12_2.drop(idx).label.values
        X12_2 = sv_data12_2.drop(idx).drop(["label"], axis=1)
        
        # Transform the data
        s12 = model1.fit_transform(X12)
        s = model1.transform(X)
        s12_2 = model1.transform(X12_2)

        # Save the model pipeline
        tf_file = open(outdir + '/cyto_sv_ml/' + str(cohort_name) + '_' + sv_type + '_' + str(k) + '_tf.pickle', 'wb')
        pickle.dump(model1, tf_file)
                    
        print(s.shape)
        print(s12.shape)
        print(s12_2.shape)
        print(np.unique(y, return_counts=True))
        print(np.unique(y12, return_counts=True))
        print(np.unique(y12_2, return_counts=True))

        print("# Train models with AutoML")
        automl = AutoML(mode="Explain", n_jobs=6
