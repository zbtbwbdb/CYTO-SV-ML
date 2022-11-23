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
import os, sys, subprocess, random, pickle, warnings, getopt

from time import time
from sklearn.utils import Bunch
from sklearn.datasets import fetch_species_distributions
import statistics
from sklearn.multiclass import OneVsRestClassifier
from itertools import cycle
from matplotlib.pyplot import *
from numpy import sqrt, argmax
from sklearn.metrics import mean_squared_error
from supervised.automl import AutoML 

warnings.filterwarnings("ignore")
#warnings.filterwarnings("ignore", category=DeprecationWarning) 

############################################################################################################################################################# 
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"s:o:n:")
for op, value in opts:
	if op == "-s":
	    sv_cohort = str(value)
	if op == "-o":
	    outdir = str(value)
	if op == "-n":
	    step_no = int(value)
sv_type_vector=['CNV', 'TRS']
subprocess.call(["sudo", "rm", "-rf", outdir])
subprocess.call(["mkdir", outdir])
kfolds=5
############################################################################################################################################################# 

for sv_type in sv_type_vector:
    # load sv data matrix with donor label
    sv_data = pd.read_csv(str(sv_cohort)+'_'+sv_type+'.csv',sep='\t', header=0, index_col=None, keep_default_na=False)
    sv_data12_2=sv_data[sv_data['label']!=0].copy()
    sv_data12_0 = sv_data[sv_data['label']==0].copy()
    sv_data12_2.index=range(sv_data12_2.shape[0])
    y_false_sv=np.array(sv_data12_2[sv_data12_2['label']== -1].index)
    y_germline_sv=np.array(sv_data12_2[sv_data12_2['label']== 1].index)
    y_somatic_sv=np.array(sv_data12_2[sv_data12_2['label']== 2].index)

    # define catagorical and continous features
    enc = preprocessing.LabelEncoder()
    categorical_varaibles=[] #['Chrom','sv_type','chr2']
    categorical_columns=[i for i in sv_data12_2.columns if i in categorical_varaibles]
    numerical_columns=[i for i in sv_data12_2.columns if i not in categorical_varaibles +['label']]
    # set up the pipeline for data transformation and modelling
    preprocessor = make_column_transformer(
         (OneHotEncoder(drop="if_binary"), categorical_columns),
        (StandardScaler(), numerical_columns),
        remainder="passthrough")
    model1 = make_pipeline(preprocessor)
                   
############################################################################################################################################################# 
                    
    # tune sub-oversampling and split training/testing/validation 
    y = sv_data12_2.label.values
    X = sv_data12_2.drop(["label"], axis=1)
    sv_ml_metrics_ts_file=open(outdir+'/'+str(sv_cohort)+'_'+sv_type+'_'+'sv_ml_metrics_sub.csv','w')  # open the text file for save model performance metrics 
    sv_ml_metrics_ts={'mic_pre':[],'mac_pre':[],'mic_rec':[],'mac_rec':[]}  
    sampling_step=round( ( (len(y_somatic_sv)*0.9+len(y_false_sv)*0.9)/3 - 2*len(y_germline_sv)*0.9/3 ) /step_no ) 
    for n in range(step_no):
        print(str(n)+": "+time.ctime())
        label_avg=round(len(y_germline_sv)*0.9)+sampling_step*n
        for m in range(kfolds):
            np.random.seed(m)
            if label_avg<=round(len(y_false_sv)*0.9):
                false_sv_idx=np.random.choice(np.random.choice(y_false_sv,round(len(y_false_sv)*0.9),replace=False),label_avg,replace=False)    
            else:
                false_sv_idx=np.random.choice(np.random.choice(y_false_sv,round(len(y_false_sv)*0.9),replace=False),label_avg,replace=True)         
            if label_avg<=round(len(y_germline_sv)*0.9):
                germline_sv_idx=np.random.choice(np.random.choice(y_germline_sv,round(len(y_germline_sv)*0.9),replace=False),label_avg,replace=False)
            else:
                germline_sv_idx=np.random.choice(np.random.choice(y_germline_sv,round(len(y_germline_sv)*0.9),replace=False),label_avg,replace=True)           
            somatic_sv_idx=np.random.choice(np.random.choice(y_somatic_sv,round(len(y_somatic_sv)*0.9),replace=False),label_avg,replace=True)
            idx=np.r_[germline_sv_idx, false_sv_idx,somatic_sv_idx]

            y12 = sv_data12_2.iloc[idx,:].label.values
            X12 = sv_data12_2.iloc[idx,:].drop(["label"], axis=1)
            y12_2 = sv_data12_2.drop(sv_data12_2.index[[idx]]).label.values
            X12_2 = sv_data12_2.drop(sv_data12_2.index[[idx]]).drop(["label"], axis=1)
            X12_0 = sv_data12_0.drop(["label"], axis=1)
            s12=model1.fit_transform(X12)
            s=model1.transform(X)
            s12_2=model1.transform(X12_2)
            s12_0=model1.transform(X12_0)

    #         print(s.shape)
    #         print(s12.shape)
    #         print(s12_2.shape)
    #         print(np.unique(y ,return_counts=True))
    #         print(np.unique(y12,return_counts=True))
    #         print(np.unique(y12_2,return_counts=True))

            # train models with AutoML
            automl = AutoML(mode="Explain", results_path=outdir+'/'+str(sv_cohort)+'_'+sv_type+'_'+str(n)+'_'+str(m)+'_ts_EXP')
            # model fitting
            automl.fit(s12, y12)

            predictions = automl.predict_all(s)
            predictions['label'].value_counts()
            print("Test accuracy:", accuracy_score(y, predictions["label"].astype(int)))
            print(predictions['label'].value_counts())
            print("Test micro precision_score:", precision_score(y, predictions["label"].astype(int),average='micro'))
            print("Test macro precision_score:", precision_score(y, predictions["label"].astype(int),average='macro'))
            print("Test micro recall_score:", recall_score(y, predictions["label"].astype(int),average='micro'))
            print("Test macro recall_score:", recall_score(y, predictions["label"].astype(int),average='macro'))

            predictions12_2 = automl.predict_all(s12_2)
            predictions12_2['label'].value_counts()
            print("Test accuracy:", accuracy_score(y12_2, predictions12_2["label"].astype(int)))
            print(predictions12_2['label'].value_counts())
            print("Test micro precision_score:", precision_score(y12_2, predictions12_2["label"].astype(int),average='micro'))
            print("Test macro precision_score:", precision_score(y12_2, predictions12_2["label"].astype(int),average='macro'))
            print("Test micro recall_score:", recall_score(y12_2, predictions12_2["label"].astype(int),average='micro'))
            print("Test macro recall_score:", recall_score(y12_2, predictions12_2["label"].astype(int),average='macro'))
            sv_ml_metrics_ts['mic_pre'].append(precision_score(y12_2, predictions12_2["label"].astype(int),average='micro'))
            sv_ml_metrics_ts['mac_pre'].append(precision_score(y12_2, predictions12_2["label"].astype(int),average='macro'))
            sv_ml_metrics_ts['mic_rec'].append(recall_score(y12_2, predictions12_2["label"].astype(int),average='micro'))
            sv_ml_metrics_ts['mac_rec'].append(recall_score(y12_2, predictions12_2["label"].astype(int),average='macro'))

            ## Display the visualization of the Confusion Matrix.
            cf_matrix = confusion_matrix(y12_2, predictions12_2['label'])
            tf=np.vectorize(lambda x: round(x,2))
            ax = sns.heatmap(tf(cf_matrix), annot=True, cmap='Blues',fmt='g')
            ax.set_title('Seaborn Confusion Matrix with labels\n\n');
            ax.set_xlabel('\nPredicted Values')
            ax.set_ylabel('Actual Values ');
            ## Ticket labels - List must be in alphabetical order
            ax.xaxis.set_ticklabels(['-1','1','2'])
            ax.yaxis.set_ticklabels(['-1','1','2'])
            ax.axhline(y=0, color='k',linewidth=1)
            ax.axhline(y=cf_matrix.shape[1], color='k',linewidth=2)
            ax.axvline(x=0, color='k',linewidth=1)
            ax.axvline(x=cf_matrix.shape[0], color='k',linewidth=2)
            sns.set(rc={'figure.figsize':(10,10)})
            plt.show()
            plt.savefig(outdir+'/'+str(sv_cohort)+'_'+sv_type+'_'+str(n)+'_'+str(m)+'_ts_model_confusion_matrix.pdf')

            # compute AUC  for model performance evaluation of multiclass
            y1 = label_binarize(y12_2, classes=[-1, 1, 2])
            y1
            fpr = dict()
            tpr = dict()
            roc_auc = dict()
            n_classes=y1.shape[1]
            lw=2
            for i in range(n_classes):
                print(i)
                fpr[i], tpr[i], _ = metrics.roc_curve(y1[:,i], predictions12_2.iloc[:,i])
                roc_auc[i] = auc(fpr[i], tpr[i])

            # Compute micro-average ROC curve and ROC area
            fpr["micro"], tpr["micro"], _ = roc_curve(y1[:,i].ravel(), predictions12_2.iloc[:,i].ravel())
            roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

            # First aggregate all false positive rates
            all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
            sv_class=['-1','1','2']

            # Then interpolate all ROC curves at this points
            mean_tpr = np.zeros_like(all_fpr)
            for i in range(n_classes):
                mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])
            # Finally average it and compute AUC
            mean_tpr /= n_classes

            fpr["macro"] = all_fpr
            tpr["macro"] = mean_tpr
            roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

            # Plot all ROC curves
            plt.figure()
            plt.plot(
                fpr["micro"],
                tpr["micro"],
                label="micro-average ROC curve (area = {0:0.2f})".format(roc_auc["micro"]),
                color="deeppink",
                linestyle=":",
                linewidth=4,
            )

            plt.plot(
                fpr["macro"],
                tpr["macro"],
                label="macro-average ROC curve (area = {0:0.2f})".format(roc_auc["macro"]),
                color="navy",
                linestyle=":",
                linewidth=4,
            )

            colors = cycle(["aqua", "darkorange", "cornflowerblue"])
            for i, color in zip(range(n_classes), colors):
                plt.plot(
                    fpr[i],
                    tpr[i],
                    color=color,
                    lw=lw,
                    label="ROC curve of class {0} (area = {1:0.2f})".format(sv_class[i], roc_auc[i]),
                )

            plt.plot([0, 1], [0, 1], "k--", lw=lw)
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel("False Positive Rate")
            plt.ylabel("True Positive Rate")
            plt.title("Some extension of Receiver operating characteristic to multiclass")
            plt.legend(loc="lower right")
            plt.show()
            plt.savefig(outdir+'/'+str(sv_cohort)+'_'+sv_type+'_'+str(n)+'_'+str(m)+'_ts_model_aucroc_curve.pdf')
    for key,value in sv_ml_metrics_ts.items():
        sv_ml_metrics_ts_file.write(str(key)+'\t'+'\t'.join(str(w) for w in value)+'\n')