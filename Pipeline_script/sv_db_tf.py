import sys,os,re
import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

# read the sv input and database vcf
SV_database_name=str(sys.argv[2])
out_vcf=str(sys.argv[1])+'_anno' 
in_vcf=pd.read_csv(out_vcf,sep='\t',header=0,keep_default_na=False)  
db_vcf=pd.read_csv(str(sys.argv[1])+'.'+str(sys.argv[2]),sep='\t',header=None,names=range(10))    

# simplify the sv database vcf
db_vcf=db_vcf.rename(columns={db_vcf.columns[6]: SV_database_name})
db_vcf.loc[db_vcf[SV_database_name]!='Not_in_database',SV_database_name]=db_vcf.iloc[:,8]
db_vcf_sim=db_vcf.iloc[:, 5:7]
#print(db_vcf_sim.iloc[0:2,:])
db_vcf_sim.columns=['sv_id',SV_database_name]
if SV_database_name not in in_vcf.columns:
    in_db_vcf=pd.merge(in_vcf, db_vcf_sim, on=['sv_id'], how='inner')
    
    # label NA database
    if  in_db_vcf[SV_database_name].isnull().any():
#    if  db_vcf_sim['db_label'].str.contains('NAN').any():        
        in_db_vcf[SV_database_name]='Not_in_database' 

# export to csv
in_db_vcf.to_csv(out_vcf,sep='\t',index=False,header=True)
