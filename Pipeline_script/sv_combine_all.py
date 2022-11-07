import sys,os,re
import pandas as pd
import numpy as np
from pathlib import Path

# read the sv input and database vcf
#info_vcf=pd.read_csv(str(sys.argv[1])+'sv.info',sep='\t',header=None,keep_default_na=False)  
#typer_vcf=pd.read_csv(str(sys.argv[1])+'svtyper.info',sep='\t',header=None,keep_default_na=False)  
db_vcf=pd.read_csv(str(sys.argv[1])+'.all_anno',sep='\t',header=None)    
sc_vcf=pd.read_csv(str(sys.argv[1])+'.bpst_bpend.fa.out.kz.index_complex',sep='\t',header=None) 
out_vcf=pd.read_csv(str(sys.argv[1])+'.all_combine',sep='\t',header=None)

# matching database label 
all_combine_vcf=pd.merge([ sc_vcf, db_vcf],on='sv_id')

# export to csv
all_combine_vcf.to_csv(out_vcf,sep='\t',index=False,header=True)
