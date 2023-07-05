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

warnings.filterwarnings("ignore")

############################################################################################################################################################# 
wd = sys.path[0]
trs_sv_cutoff=1000 # trs_sv breakpoint distance based cutoff
nontrs_sv_cutoff=0.9 # nontrs_sv sequence overlapping based cutoff
opts,args = getopt.getopt(sys.argv[1:],"s:o:t:c:")
for op, value in opts:
	if op == "-s":
	    cohort_name = str(value)
	if op == "-o":
	    outdir = str(value)
	if op == "-t":
	    trs_sv_cutoff = int(value)
	if op == "-c":
	    nontrs_sv_cutoff = float(value)  
sv_type_vector=['trs','nontrs']

############################################################################################################################################################# 
sv_data = pd.read_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all',sep='\t', header=0, index_col=None, keep_default_na=False)
print(sv_data.columns)
for sv_type in sv_type_vector:
    if sv_type=='trs':
        sv_data_01=sv_data[~(sv_data['sv_type'].isin(['DEL','DUP','INV','INS']))].copy()  
        sv_data_01.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_trs_o',sep='\t', header=True, index=None) 
        sv_data_tf=sv_dataframe_transform.trs_sv_data_transform(sv_data_01,trs_sv_cutoff)
        sv_data_1=sv_data_tf[0]
        sv_data_1.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_trs',sep='\t', header=True, index=None)
        sv_data_10=sv_data_tf[1]
        sv_data_10.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_trs_all',sep='\t', header=True, index=None)        
    else:
        sv_data_01=sv_data[sv_data['sv_type'].isin(['DEL','DUP','INV','INS'])].copy()     
        sv_data_01.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_nontrs_o',sep='\t', header=True, index=None) 
        sv_data_tf=sv_dataframe_transform.nontrs_sv_data_transform(sv_data_01,nontrs_sv_cutoff)        
        sv_data_1=sv_data_tf[0]
        sv_data_1.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_nontrs',sep='\t', header=True, index=None)        
        sv_data_10=sv_data_tf[1]
        sv_data_10.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_nontrs_all',sep='\t', header=True, index=None)            
    sv_summary_plot=sv_dataframe_transform.sv_data_summary_plot(sv_data_1)    
    sv_summary_plot.savefig(outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+'sv_summary_plot.pdf', transparent=True)  # open the plot file for save sv data summary   
