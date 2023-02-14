import sys,getopt,os
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
from sklearn.ensemble import GradientBoostingRegressor,GradientBoostingClassifier, AdaBoostClassifier
from sklearn.metrics import mean_squared_error, median_absolute_error, roc_curve, auc, accuracy_score, precision_score, recall_score, precision_recall_curve, confusion_matrix
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer, TransformedTargetRegressor
from sklearn.pipeline import make_pipeline
from supervised.automl import AutoML 

#parameter setting
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"t:n:i:")
inFile = ""

for op, value in opts:
	if op == "-t":
	    trs_dir = str(value)
	if op == "-n":
	    nontrs_dir = str(value)
	if op == "-i":
	    inFile =  str(value)

if inFile == "":
	print("-i invalid")
	sys.exit()
      
# load the optimal model
automl_trs = AutoML(mode="Explain", n_jobs= 6, results_path=trs_dir+'_EXP')
automl_nontrs = AutoML(mode="Explain", n_jobs= 6, results_path=nontrs_dir+'_EXP')

# load the training data tranformer
model_trs=pickle.load(open(trs_dir+'.pickle', 'rb'))
model_nontrs=pickle.load(open(nontrs_dir+'.pickle', 'rb')) 

# read in all WGS SV data for model prediction
X_trs=pd.read_csv(inFile+'.sv.all.combine_all_trs',sep='\t',header=0, index_col=None, keep_default_na=False)
X_nontrs=pd.read_csv(inFile+'.sv.all.combine_all_nontrs',sep='\t',header=0, index_col=None, keep_default_na=False)
X_all=pd.read_csv(inFile+'.sv.all.combine_all',sep='\t',header=0, index_col=None, keep_default_na=False)

# transform all the original SV data
s_trs=model_trs.transform(X_trs)
s_nontrs=model_nontrs.transform(X_nontrs)

# predict all the transformed SV data
trs_label=automl_trs.predict_all(s_trs)
nontrs_label=automl_nontrs.predict_all(s_nontrs)

# attach the prediction label to all the original SV data
X_trs['predict_label','prediction_TA','prediction_TG','prediction_TS']=trs_label
X_nontrs['predict_label','prediction_TA','prediction_TG','prediction_TS']=nontrs_label
tmp_id=['sv_type','sv_bp_st_cc1','sv_bp_end_cc1','sv_bp_st_cc_v1','sv_bp_end_cc_v1','sv_bp_st_cc_v2','sv_bp_end_cc_v2','cipos_range','ciend_range','PR_read_ratio','SR_read_ratio']
X_trs['tmp_id']=X_trs.loc[:,tmp_id].apply(lambda x: ':'.join(str(w) for w in x), axis=1)      
X_nontrs['tmp_id']=X_nontrs.loc[:,tmp_id].apply(lambda x: ':'.join(str(w) for w in x), axis=1) 
X_tmp=pd.concat([X_trs['tmp_id','predict_label','prediction_TA','prediction_TG','prediction_TS'],X_nontrs['tmp_id','predict_label','prediction_TA','prediction_TG','prediction_TS']])
X_all['tmp_id']=X_all.loc[:,tmp_id].apply(lambda x: ':'.join(str(w) for w in x), axis=1) 
X_all['predict_label','prediction_TA','prediction_TG','prediction_TS']=X_tmp['predict_label','prediction_TA','prediction_TG','prediction_TS']

# re-organize the data
X_all=X_all.rename(columns={'old_col':'new_col','old_col':'new_col'})
col_drop=[str(w) for w in X_trs.columns if w not in ['sv_type', 'sv_chr', 'sv_chr2', 'sv_read_r', 'sv_read_a', 'sv_read_ratio', 'sv_read_diff', 'sv_bp_st', 'sv_bp_end', 'sv_bp_st_ci0', 'sv_bp_st_ci1', 'sv_bp_end_ci0', 'sv_bp_end_ci1', 'sv_bp_st_ci_range', 'sv_bp_end_ci_range', 'sv_bp_st_cc_v1', 'sv_bp_end_cc_v1', 'sv_database', 'predict_label', 'prediction_TA', 'prediction_TG', 'prediction_TS', 'label', 'predict_max']]
X_all=X_all.drop(col_drop,axis=1)

# save all the SV data
pd.to_csv(X_trs,inFile+'trs',sep='\t',index=False,header=True)
pd.to_csv(X_nontrs,inFile+'nontrs',sep='\t',index=False,header=True)
