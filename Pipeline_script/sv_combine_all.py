import sys,os,re
import pandas as pd
import numpy as np
from pathlib import Path

# read the sv input and database vcf
SV_database_name=str(sys.argv[2])
info_vcf=pd.read_csv(str(sys.argv[1])+'svtyper.info',sep='\t',header=True,keep_default_na=False)  
sc_vcf=pd.read_csv(str(sys.argv[1])+'.all_anno',sep='\t',header=None)    
db_vcf=pd.read_csv(str(sys.argv[1])+'.bpst_bpend.fa.out.kz.index_complex',sep='\t',header=None) 
out_vcf=pd.read_csv(str(sys.argv[1])+'.all_combine',sep='\t',header=None)

# matching database label 
all_combine_vcf=pd.merge([info_vcf, sc_vcf, db_vcf],on='sv_id')

# export to csv
all_combine_vcf.to_csv(out_vcf,sep='\t',index=False,header=True)
