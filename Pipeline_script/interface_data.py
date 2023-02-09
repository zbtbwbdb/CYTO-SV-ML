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
X_trs['label']=trs_label['label']
X_nontrs['label']=nontrs_label['label']

# save all the SV data
pd.to_csv(X_trs,inFile+'trs',sep='\t',index=False,header=True)
pd.to_csv(X_nontrs,inFile+'nontrs',sep='\t',index=False,header=True)
