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
opts,args = getopt.getopt(sys.argv[1:],"s:o:k:")
for op, value in opts:
	if op == "-s":
	    cohort_name = str(value)
	if op == "-o":
	    outdir = str(value)
	if op == "-k":
	    kfolds = int(value)
sv_type_vector=['trs','nontrs']
# subprocess.call(["sudo","mkdir", "-p", outdir,'/cyto_sv_ml'])

############################################################################################################################################################# 
sv_data = pd.read_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all',sep='\t', header=0, index_col=None, keep_default_na=False)
print(sv_data.columns)
for sv_type in sv_type_vector:
    if sv_type=='trs':
        sv_data_01=sv_data[~(sv_data['sv_type'].isin(['DEL','DUP','INV','INS']))].copy()  
        sv_data_01.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_trs_o',sep='\t', header=True, index=None)        
        sv_data_1=sv_dataframe_transform.trs_sv_data_transform(sv_data_01)[0]
        sv_data_1.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_trs',sep='\t', header=True, index=None)
        sv_data_10=sv_dataframe_transform.trs_sv_data_transform(sv_data_01)[1]
        sv_data_10.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_trs_all',sep='\t', header=True, index=None)        
    else:
        sv_data_01=sv_data[sv_data['sv_type'].isin(['DEL','DUP','INV','INS'])].copy()     
        sv_data_01.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_nontrs_o',sep='\t', header=True, index=None)       
        sv_data_1=sv_dataframe_transform.nontrs_sv_data_transform(sv_data_01)[0]
        sv_data_1.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_nontrs',sep='\t', header=True, index=None)        
        sv_data_10=sv_dataframe_transform.nontrs_sv_data_transform(sv_data_01)[1]
        sv_data_10.to_csv(outdir+'/'+str(cohort_name)+'.sv.all.combine_all_nontrs_all',sep='\t', header=True, index=None)            
    sv_summary_plot=sv_dataframe_transform.sv_data_summary_plot(sv_data_1)    
    sv_summary_plot.savefig(outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+'sv_summary_plot.pdf', transparent=True)  # open the plot file for save sv data summary   

    print("# load sv data matrix with label")
    sv_data12_2=sv_data_1[sv_data_1['label']!=0].copy()
    sv_data12_2=sv_data12_2.drop(['sv_type'],axis=1)
    sv_data12_0 = sv_data_1[sv_data_1['label']==0].copy()
    sv_data12_2.index=range(sv_data12_2.shape[0])
    y_false_sv=np.array(sv_data12_2[sv_data12_2['label']== -1].index)
    y_germline_sv=np.array(sv_data12_2[sv_data12_2['label']== 1].index)
    y_somatic_sv=np.array(sv_data12_2[sv_data12_2['label']== 2].index)

    print("# define catagorical and continous features")
    enc = preprocessing.LabelEncoder()
    categorical_varaibles=[] #['Chrom','sv_type','chr2']
    categorical_columns=[i for i in sv_data12_2.columns if i in categorical_varaibles]
    numerical_columns=[i for i in sv_data12_2.columns if i not in categorical_varaibles +['label']]
    print("# set up the pipeline for data transformation and modelling")
    preprocessor = make_column_transformer(
         (OneHotEncoder(drop="if_binary"), categorical_columns),
        (StandardScaler(), numerical_columns),
        remainder="passthrough")
    model1 = make_pipeline(preprocessor)
                   
############################################################################################################################################################# 
                    
    print("# tune sub-oversampling and split training/testing/validation") 
    y = sv_data12_2.label.values
    X = sv_data12_2.drop(["label"], axis=1)
    sv_ml_metrics_ts_file=open(outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+'sv_ml_metrics_sub.csv','w')  # open the text file for save model performance metrics 
    sv_ml_metrics_ts={'mic_pre':[],'mac_pre':[],'mic_rec':[],'mac_rec':[]}  
    label_avg=median([round(len(y_somatic_sv)*0.9),round(len(y_germline_sv)*0.9),round(len(y_false_sv)*0.9)])
    for k in range(kfolds):
        np.random.seed(k*2)
        if label_avg<=round(len(y_false_sv)*0.9):
            false_sv_idx=np.random.choice(np.random.choice(y_false_sv,round(len(y_false_sv)*0.9),replace=False),label_avg,replace=False)    
        else:
            false_sv_idx=np.random.choice(np.random.choice(y_false_sv,round(len(y_false_sv)*0.9),replace=False),label_avg,replace=True)         
        if label_avg<=round(len(y_germline_sv)*0.9):
            germline_sv_idx=np.random.choice(np.random.choice(y_germline_sv,round(len(y_germline_sv)*0.9),replace=False),label_avg,replace=False)
        else:
            germline_sv_idx=np.random.choice(np.random.choice(y_germline_sv,round(len(y_germline_sv)*0.9),replace=False),label_avg,replace=True)                   
        if label_avg<=round(len(y_somatic_sv)*0.9):
            somatic_sv_idx=np.random.choice(np.random.choice(y_somatic_sv,round(len(y_somatic_sv)*0.9),replace=False),label_avg,replace=False)
        else:                
            somatic_sv_idx=np.random.choice(np.random.choice(y_somatic_sv,round(len(y_somatic_sv)*0.9),replace=False),label_avg,replace=True)
        idx=np.r_[germline_sv_idx, false_sv_idx,somatic_sv_idx]

        y12 = sv_data12_2.iloc[idx,:].label.values
        X12 = sv_data12_2.iloc[idx,:].drop(["label"], axis=1)
        y12_2 = sv_data12_2.drop(idx).label.values
        X12_2 = sv_data12_2.drop(idx).drop(["label"], axis=1)
        #X12_0 = sv_data12_0.drop(["label"], axis=1)
        s12=model1.fit_transform(X12)
        s=model1.transform(X)
        s12_2=model1.transform(X12_2)
        #s12_0=model1.transform(X12_0)
        # open a file, where you ant to store the data
        tf_file = open(outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_tf.pickle', 'wb')
        pickle.dump(model1,tf_file)
                    
        print(s.shape)
        print(s12.shape)
        print(s12_2.shape)
        print(np.unique(y ,return_counts=True))
        print(np.unique(y12,return_counts=True))
        print(np.unique(y12_2,return_counts=True))

        print("# train models with AutoML")
        automl = AutoML(mode="Explain", n_jobs= 6, results_path=outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_ts_EXP')
        # model fitting
        automl.fit(s12, y12)
        
        print(" #copy shap feature importance plot")
        shap_plot=outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_ts_EXP/'+'*_Default_Xgboost/learner_fold_0_shap_summary.png'
        shap_plot2=outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_ts_EXP/learner_fold_0_shap_summary.png'
        os. system("sudo cp "+ shap_plot + " " + shap_plot2) 
        
        print("# model performance summary metrics") 
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

        print("## Display the visualization of the Confusion Matrix.")
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
        plt.savefig(outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_ts_model_confusion_matrix.pdf')

        print("# compute AUC  for model performance evaluation of multiclass")
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

        print("# Plot all ROC curves")
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
        plt.savefig(outdir+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_ts_model_aucroc_curve.pdf')
    for key,value in sv_ml_metrics_ts.items():
        sv_ml_metrics_ts_file.write(str(key)+','+','.join(str(w) for w in value)+'\n')
