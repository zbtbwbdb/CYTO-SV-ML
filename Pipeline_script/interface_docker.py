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

# load the training data tranformer

# read in all WGS SV data

# transform all the original SV data

# predict all the transformed SV data

# attach the prediction label to all the original SV data

# save all the SV data

