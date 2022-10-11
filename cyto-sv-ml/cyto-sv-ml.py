import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import seaborn as sns
from pandas import read_csv
from numpy import asarray
from itertools import cycle
from numpy import sqrt, argmax
from itertools import *
import random
import pickle
from time import time
import statistics
from sklearn.utils import Bunch
from supervised.automl import AutoML # mljar-supervised

from sklearn import metrics, datasets
from sklearn.metrics import mean_squared_error, median_absolute_error, roc_curve, auc, accuracy_score, precision_score, recall_score, precision_recall_curve, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.pipeline import make_pipeline
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler, OneHotEncoder, OrdinalEncoder, label_binarize, LabelEncoder
from sklearn.ensemble import GradientBoostingRegressor, GradientBoostingClassifier, AdaBoostClassifier


# Input data
sv=pd.read_csv("input_sv.vcf",sep="\t", header=0, index_col=None, keep_default_na=False)
sv_type=sv['sv_type']

if (sv_type=='BND'):
    cyto_sv_ml= AutoML(mode="Explain", algorithms=['Xgboost'],results_path="trs-cyto-sv")
    sv_columns=pd.read_csv("trs_col.index",sep="\t", header=None, index_col=None, keep_default_na=False)
    data_tf=pickle.load(open('trs_tf.pickle', 'rb'))    
else:
    cyto_sv_ml= AutoML(mode="Explain", algorithms=['Xgboost'],results_path="cnv-cyto-sv")
    sv_columns=pd.read_csv("cnv_col.index",sep="\t", header=None, index_col=None, keep_default_na=False)
    data_tf=pickle.load(open('cnv_tf.pickle', 'rb'))

# data transformation setting for data_tf
# preprocessor = make_column_transformer((OneHotEncoder(drop="if_binary"), categorical_columns),(StandardScaler(), numerical_columns),remainder="passthrough")

# data transformation
sv_sel=sv.iloc[:,np.r_[np.where(sv.columns.isin(sv_columns.iloc[0,:]))]
X = sv_sel.drop("label", axis=1)               
y = sv_sel.label.values
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,stratify=X[['label']], random_state=0)
X_train_tf=data_tf.transform(X_train)  
X_test_tf=data_tf.transform(X_test)                 
               
# AUTOML model fit with training data and validate with testing data              
sv_train_model = cyto_sv_ml.fit(X_train_tf)
predictions_test = sv_train_model.predict_all(X_test_tf)
               
# model performance metrics
print("Test accuracy:", accuracy_score(y_test, predictions_test["label"].astype(int)))
print(predictions_test['label'].value_counts())
print("Test micro precision_score:", precision_score(y_test, predictions_test["label"].astype(int),average='micro'))
print("Test macro precision_score:", precision_score(y_test, predictions_test["label"].astype(int),average='macro'))
print("Test micro recall_score:", recall_score(y_testing, predictions_testing["label"].astype(int),average='micro'))
print("Test macro recall_score:", recall_score(y_testing, predictions_testing["label"].astype(int),average='macro'))
