import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import read_csv
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder
from numpy import asarray
from sklearn import svm, metrics, datasets
from sklearn.metrics import median_absolute_error, roc_curve, auc, accuracy_score, precision_score, recall_score, precision_recall_curve, confusion_matrix
from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import OneHotEncoder, label_binarize
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import GradientBoostingClassifier, AdaBoostClassifier
from sklearn.preprocessing import StandardScaler
from itertools import *
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import Ridge
from sklearn.compose import TransformedTargetRegressor
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder
from sklearn import preprocessing
import random
import pickle
from time import time
from sklearn.utils import Bunch
from sklearn.datasets import fetch_species_distributions
import statistics
from sklearn.multiclass import OneVsRestClassifier
from itertools import cycle
from matplotlib.pyplot import *
from numpy import sqrt, argmax
from sklearn.metrics import mean_squared_error
from supervised.automl import AutoML # mljar-supervised

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
sv_tf=data_tf.transform(sv)
sv_train_model = cyto_sv_ml.fit(sv_tf)
print(sv_train_model)
               
# model performance metrics
print("Test accuracy:", accuracy_score(y_test, predictions_test["label"].astype(int)))
print(predictions_test['label'].value_counts())
print("Test micro precision_score:", precision_score(y_test, predictions_test["label"].astype(int),average='micro'))
print("Test macro precision_score:", precision_score(y_test, predictions_test["label"].astype(int),average='macro'))
print("Test micro recall_score:", recall_score(y_test, predictions_test["label"].astype(int),average='micro'))
print("Test macro recall_score:", recall_score(y_test, predictions_test["label"].astype(int),average='macro'))
