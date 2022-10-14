import os, sys, statistics, random, pickle
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
from time import time
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
sv=pd.read_csv(sys.argv[1],sep="\t", header=0, index_col=None, keep_default_na=False)
sv_type=sv['sv_type']

if (sv_type=='BND'):
    cyto_sv_ml= AutoML(mode="Explain", results_path="trs-cyto-sv")
    sv_columns=pd.read_csv("trs_col.index",sep="\t", header=None, index_col=None, keep_default_na=False)
    data_tf=pickle.load(open('trs_tf.pickle', 'rb'))    
else:
    cyto_sv_ml= AutoML(mode="Explain", results_path="cnv-cyto-sv")
    sv_columns=pd.read_csv("cnv_col.index",sep="\t", header=None, index_col=None, keep_default_na=False)
    data_tf=pickle.load(open('cnv_tf.pickle', 'rb'))

# data transformation setting for data_tf
# categorical_columns=[i for i in sv_columns if i in categorical_varaibles]
# numerical_columns=[i for i in sv_columns if i not in categorical_varaibles +['label']]
# preprocessor = make_column_transformer((OneHotEncoder(drop="if_binary"), categorical_columns),(StandardScaler(), numerical_columns),remainder="passthrough")
# data_tf = make_pipeline(preprocessor)
# Class -1 is for systematic artifact SVs; Class 1 is for true germline SVs; Class 2 is for true cytogenetic somatic SVs

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
               
# Model performance metrics
print("Test accuracy:", accuracy_score(y_test, predictions_test["label"].astype(int)))
print(predictions_test['label'].value_counts())
print("Test micro precision_score:", precision_score(y_test, predictions_test["label"].astype(int),average='micro'))
print("Test macro precision_score:", precision_score(y_test, predictions_test["label"].astype(int),average='macro'))
print("Test micro recall_score:", recall_score(y_testing, predictions_testing["label"].astype(int),average='micro'))
print("Test macro recall_score:", recall_score(y_testing, predictions_testing["label"].astype(int),average='macro'))

## Display the confusion matrix
cf_matrix = confusion_matrix(y_test, predictions_test['label'],normalize="true")
tf=np.vectorize(lambda x: round(x,2))
ax = sns.heatmap(tf(cf_matrix), annot=True, cmap='Blues',fmt='g')               
ax.set_title('Seaborn Confusion Matrix with labels\n\n');
ax.set_xlabel('\nPredicted Values')
ax.set_ylabel('Actual Values ');
ax.xaxis.set_ticklabels(['-1','1','2'])
ax.yaxis.set_ticklabels(['-1','1','2'])
ax.axhline(y=0, color='k',linewidth=1)
ax.axhline(y=cf_matrix.shape[1], color='k',linewidth=2)
ax.axvline(x=0, color='k',linewidth=1)
ax.axvline(x=cf_matrix.shape[0], color='k',linewidth=2)
sns.set(rc={'figure.figsize':(10,10)})
plt.show()       
           
# AUCROC curve  for model performance evaluation
y1 = label_binarize(y_test, classes=[-1, 1, 2])
fpr = dict()
tpr = dict()
roc_auc = dict()
n_classes=y1.shape[1]
lw=2
for i in range(n_classes):
    print(i)
    fpr[i], tpr[i], _ = metrics.roc_curve(y1[:,i], predictions_test.iloc[:,i])
    roc_auc[i] = auc(fpr[i], tpr[i])
colors = cycle(['blue', 'red', 'green'])
for i, color in zip(range(n_classes), colors):
    plt.plot(fpr[i], tpr[i], color=color, lw=2,
             label='ROC curve of class {0} (area = {1:0.2f})'
             ''.format(i, roc_auc[i]))
plt.plot([0, 1], [0, 1], 'k--', lw=lw)
plt.xlim([-0.05, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic for multi-class data')
plt.legend(loc="lower right")
plt.show()  
