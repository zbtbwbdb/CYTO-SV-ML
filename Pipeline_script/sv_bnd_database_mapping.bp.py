################################################################################################################
# TRS SV simplified data format
# ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.nobnd
# ---example-----------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_chr2  sv_type  sv_id
# chr1    100000       1000000    chr21    BND      chr1:100000:1000000:chr21:BND:MantaBND***** 
# ----------------------------------------------------------------------------
# TRS SV simplified database format
# ${main_dir}/SV_database/${SV_database_name}.gz
# ---example------------------------------------------------------------------
# sv_chr1  sv_start_bp  sv_end_bp  sv_chr2  sv_type  database_AF
# chr1     100000       1000000    chr21    BND      0.05
# ----------------------------------------------------------------------------
#################################################################################################################

import sys,os,re
import pandas as pd

ref_data=pd.read_csv(sys.argv[1],sep='\t',header=None)#,skiprows=1)
in_vcf=pd.read_csv(sys.argv[2],sep='\t',header=0,keep_default_na=False)
out_vcf=str(sys.argv[2])+"."+str(sys.argv[3])
bp_dis=int(str(sys.argv[3]).split('_')[2])

in_vcf['database']='Not_in_database'
for i in range(in_vcf.shape[0]):
    chr1=str(in_vcf.iloc[i,0])
    bp1=int(in_vcf.iloc[i,1])
    chr2=str(in_vcf.iloc[i,3])   
    bp2=int(in_vcf.iloc[i,2])  
#     if i%1000==0:
#         print(str(i)+":\t"+chr1+"\t"+str(bp1)+"\t"+chr2+"\t"+str(bp2))
    bnd_dict=ref_data[(ref_data.iloc[:,0]==chr1) & (ref_data.iloc[:,3]==chr2) & (ref_data.iloc[:,1].astype(int)<=bp1+bp_dis)  & (ref_data.iloc[:,2].astype(int)>=(bp1-bp_dis)) & (ref_data.iloc[:,4].astype(int)<=bp2+bp_dis) & (ref_data.iloc[:,5].astype(int)>=(bp2-bp_dis))] 
    if bnd_dict.empty or bnd_dict.shape[0]==0:
        continue
    else:
        bnd_dict.iloc[:,1]=bnd_dict.iloc[:,1].astype(int)-bp1
        bnd_dict.iloc[:,2]=bp1-bnd_dict.iloc[:,2].astype(int) 
        bnd_dict.iloc[:,3]=0    
        bnd_dict.iloc[:,4]=bnd_dict.iloc[:,4].astype(int)-bp1
        bnd_dict.iloc[:,5]=bp1-bnd_dict.iloc[:,5].astype(int) 
        bnd_dict['info']=bnd_dict.apply(lambda x: max(x[1:6]),axis=1)  
        in_vcf.iloc[i,in_vcf.shape[1]-1]=min(bnd_dict['info']) 
in_vcf.to_csv(out_vcf,sep='\t',index=False,header=False)
