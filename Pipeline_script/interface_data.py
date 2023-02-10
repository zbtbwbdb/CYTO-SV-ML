import sys,getopt,os,commands,copy,subprocess,SVCNV,SVCNV_sim,SVCNV_set
#parameter setting
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"i:d:p:t:o:")
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
automl_trs = AutoML(mode="Explain", n_jobs= 6, results_path=trs_dir)
automl_nontrs = AutoML(mode="Explain", n_jobs= 6, results_path=nontrs_dir)

# load the training data tranformer
model_trs=pickle.load(open(trs_dir+'trs.pickle', 'rb'))
model_nontrs=pickle.load(open(nontrs_dir+'nontrs.pickle', 'rb'))	

# read in all WGS SV data
X_trs=pd.read_csv(inFile+'trs',sep='\t',header=0, index_col=None, keep_default_na=False)
X_nontrs=pd.read_csv(inFile+'nontrs',sep='\t',header=0, index_col=None, keep_default_na=False)

# transform all the original SV data
s_trs=model_trs.transform(X_trs)
s_nontrs=model_nontrs.transform(X_nontrs)

# predict all the transformed SV data
trs_label=automl_trs.predict_all(s_trs)
nontrs_label=automl_nontrs.predict_all(s_nontrs)

# attach the prediction label to all the original SV data
X_trs['predict_label','prediction_TA''prediction_TG''prediction_TS']=trs_label
X_nontrs['predict_label','prediction_TA''prediction_TG''prediction_TS']=nontrs_label

# re-organize the data
X_trs=X_trs.rename(columns={'old_col':'new_col','old_col':'new_col'})
col_drop=[str(w) for w in X_trs.columns if w not in ['sv_type', 'sv_chr', 'sv_chr2', 'sv_read_r', 'sv_read_a', 'sv_read_ratio', 'sv_read_diff', 'sv_bp_st', 'sv_bp_end', 'sv_bp_st_ci0', 'sv_bp_st_ci1', 'sv_bp_end_ci0', 'sv_bp_end_ci1', 'sv_bp_st_ci_range', 'sv_bp_end_ci_range', 'sv_bp_st_cc_v1', 'sv_bp_end_cc_v1', 'sv_database', 'predict_label', 'prediction_TA', 'prediction_TG', 'prediction_TS', 'label', 'predict_max']]
X_trs=X_trs.drop(col_drop,axis=1)
X_nontrs=X_nontrs.rename(columns={'old_col':'new_col','old_col':'new_col'})
col_drop=[str(w) for w in X_nontrs.columns if w not in ['sv_type', 'sv_chr', 'sv_chr2', 'sv_read_r', 'sv_read_a', 'sv_read_ratio', 'sv_read_diff', 'sv_bp_st', 'sv_bp_end', 'sv_bp_st_ci0', 'sv_bp_st_ci1', 'sv_bp_end_ci0', 'sv_bp_end_ci1', 'sv_bp_st_ci_range', 'sv_bp_end_ci_range', 'sv_bp_st_cc_v1', 'sv_bp_end_cc_v1', 'sv_database', 'predict_label', 'prediction_TA', 'prediction_TG', 'prediction_TS', 'label', 'predict_max']]
X_nontrs=X_nontrs.drop(col_drop,axis=1)

# save all the SV data
pd.to_csv(X_trs,inFile+'trs',sep='\t',index=False,header=True)
pd.to_csv(X_nontrs,inFile+'nontrs',sep='\t',index=False,header=True)
