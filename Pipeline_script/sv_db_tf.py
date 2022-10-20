import sys,os,re
import pandas as pd
from pathlib import Path

# read the sv input and database vcf
SV_database_name=str(sys.argv[2])
in_vcf=pd.read_csv(str(sys.argv[1]),sep='\t',header=True,keep_default_na=False)  
db_vcf=pd.read_csv(str(sys.argv[1])+'_'+str(sys.argv[2]),sep='\t',header=None)    
out_vcf=str(sys.argv[1])+'.anno' 

# simplify the sv database vcf
db_vcf_sim=db_vcf.iloc[:, 1:(db_vcf.shape[1]-1)]
db_vcf_sim.columns=['sv_id','db_label']
db_vcf_sim=db_vcf_sim[db_vcf_sim['db_label']=='FAIL']

# matching database label 
in_vcf[SV_database_name]='Not_in_database'
in_vcf[in_vcf['sv_id'].isin(db_vcf_sim['sv_id']),SV_database_name]='In_database'

# export to csv
in_vcf.to_csv(out_vcf,sep='\t',index=False,header=True)
