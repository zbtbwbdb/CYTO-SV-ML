import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random 
from itertools import repeat
import os, sys, getopt

############################################################################################################################################################# 

def trs_sv_data_transform(sv_data_trs,trs_sv_cutoff):
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_trs.shape[0])), k=sv_data_trs.shape[0]) 
    sv_data_trs=sv_data_trs.reset_index(drop=True)
    sv_data_trs=sv_data_trs.reindex(new_index)
    sv_data_trs=sv_data_trs.reset_index(drop=True)   
    print(sv_data_trs.shape) 
    
    # transform database label
    sv_db_list=['1000_gall','1000_g','gnomad_gall','gnomad_g','control_gall', 'control_g', 'cytoatlas_s','cosmic_s','gnomad_qc', 'dgv_g', 'centromere_qc','uwstl_s' ]
    for sv_db in sv_db_list:
#        sv_data_trs.loc[np.r_[np.where(sv_data_trs[sv_db].isnull()) ],sv_db]="Not_in_database"
        sv_data_trs.loc[np.r_[np.where(sv_data_trs[sv_db].isnull()) ],sv_db]=trs_sv_cutoff*10000
        sv_data_trs.loc[np.r_[np.where(sv_data_trs[sv_db]=='NAN')],sv_db]=trs_sv_cutoff*10000
        sv_data_trs.loc[np.r_[np.where(sv_data_trs[sv_db]=='Not_in_database')],sv_db]=trs_sv_cutoff*10000   
        
    sv_data_trs.loc[sv_data_trs['1000_gall'].astype(int)<=trs_sv_cutoff,'1000_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['1000_g'].astype(int)<=trs_sv_cutoff,'1000_g']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_gall'].astype(int)<=trs_sv_cutoff,'gnomad_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_g'].astype(int)<=trs_sv_cutoff,'gnomad_g']="YVID"
    sv_data_trs.loc[sv_data_trs['control_gall'].astype(int)<=trs_sv_cutoff,'control_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['control_g'].astype(int)<=trs_sv_cutoff,'control_g']="YVID"
    sv_data_trs.loc[sv_data_trs['cytoatlas_s'].astype(int)<=trs_sv_cutoff,'cytoatlas_s']="YVID"
    sv_data_trs.loc[sv_data_trs['cosmic_s'].astype(int)<=trs_sv_cutoff,'cosmic_s']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_qc'].astype(int)<=trs_sv_cutoff,'gnomad_qc']="YVID"
    sv_data_trs.loc[sv_data_trs['dgv_g'].astype(int)<=trs_sv_cutoff,'dgv_g']="YVID"
    sv_data_trs.loc[sv_data_trs['centromere_qc'].astype(int)<=trs_sv_cutoff,'centromere_qc']="YVID"
    sv_data_trs.loc[sv_data_trs['uwstl_s'].astype(int)<=trs_sv_cutoff,'uwstl_s']="YVID"
#     print(sv_data_trs['1000_gall'].value_counts())
#     print(sv_data_trs['1000_g'].value_counts())
#     print(sv_data_trs['gnomad_gall'].value_counts())
#     print(sv_data_trs['gnomad_g'].value_counts())
#     print(sv_data_trs['control_gall'].value_counts())
#     print(sv_data_trs['control_g'].value_counts())
#     print(sv_data_trs['cytoatlas_s'].value_counts())
#     print(sv_data_trs['cosmic_s'].value_counts())
#     print(sv_data_trs['gnomad_qc'].value_counts())
#     print(sv_data_trs['centromere_qc'].value_counts())
#     print(sv_data_trs['dgv_g'].value_counts())
#     print(sv_data_trs['uwstl_s'].value_counts())
    
    # transform sv_info variables
    sv_data_trs['ci_pos0']=0
    sv_data_trs['ci_pos1']=0
    sv_data_trs['ci_end']=0
    sv_data_trs['PR_ref']=0
    sv_data_trs['PR_alt']=0
    sv_data_trs['SR_ref']=0
    sv_data_trs['SR_alt']=0
    sv_data_trs['ci_pos0']=sv_data_trs['CIPOS'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['ci_pos1']=sv_data_trs['CIPOS'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_trs['sv_bp_end_POS']=sv_data_trs['CIEND'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['sv_bp_end_END']=sv_data_trs['CIEND'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0) 

    # read statistics transformation
    sv_data_trs.loc[(sv_data_trs['RP'].isnull()) | (sv_data_trs['RP']=='NA') | (sv_data_trs['RP']=='.'),'RP']=0  
    sv_data_trs.loc[(sv_data_trs['AP'].isnull()) | (sv_data_trs['AP']=='NA') | (sv_data_trs['AP']=='.'),'AP']=0  
    sv_data_trs.loc[(sv_data_trs['RS'].isnull()) | (sv_data_trs['RS']=='NA') | (sv_data_trs['RS']=='.'),'RS']=0          
    sv_data_trs.loc[(sv_data_trs['AS'].isnull()) | (sv_data_trs['AS']=='NA') | (sv_data_trs['AS']=='.'),'AS']=0    
    sv_data_trs.loc[(sv_data_trs['DR'].isnull()) | (sv_data_trs['DR']=='NA') | (sv_data_trs['DR']=='.'),'DR']=0  
    sv_data_trs.loc[(sv_data_trs['DV'].isnull()) | (sv_data_trs['DV']=='NA') | (sv_data_trs['DV']=='.'),'DV']=0  
    sv_data_trs.loc[(sv_data_trs['RR'].isnull()) | (sv_data_trs['RR']=='NA') | (sv_data_trs['RR']=='.'),'RR']=0          
    sv_data_trs.loc[(sv_data_trs['RV'].isnull()) | (sv_data_trs['RV']=='NA') | (sv_data_trs['RV']=='.'),'RV']=0  
    sv_data_trs.to_csv("/home/tzhang/test0",sep="\t",header=True,index=False)     
#     print(sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['AP']!=0),['RP','AP']].shape)
    sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['AP']!=0),'PR']=sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['AP']!=0),['RP','AP']].apply(lambda x: str(x['RP'])+','+str(x['AP']), axis=1)
    sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(',')) | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['AS']!=0),'SR']=sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(',')) | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['AS']!=0),['RS','AS']].apply(lambda x: str(x['RS'])+','+str(x['AS']), axis=1)             
    sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['DV']!=0),'PR']=sv_data_trs.loc[((~sv_data_trs['PR'].str.contains(',')) | (sv_data_trs['PR']=='0,0')) & (sv_data_trs['DV']!=0),['DR','DV']].apply(lambda x: str(x['DR'])+','+str(x['DV']), axis=1)
    sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(','))  | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['RV']!=0),'SR']=sv_data_trs.loc[((~sv_data_trs['SR'].str.contains(','))  | (sv_data_trs['SR']=='0,0')) & (sv_data_trs['RV']!=0),['RR','RV']].apply(lambda x: str(x['RR'])+','+str(x['RV']), axis=1)  
    sv_data_trs.loc[~sv_data_trs['PR'].str.contains(','),'PR']='0,0'
    sv_data_trs.loc[~sv_data_trs['SR'].str.contains(','),'SR']='0,0'         
         
    sv_data_trs['PR_ref']=sv_data_trs['PR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['PR_alt']=sv_data_trs['PR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_trs['SR_ref']=sv_data_trs['SR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_trs['SR_alt']=sv_data_trs['SR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)        

    # replace na
    sv_data_trs=sv_data_trs.fillna(0)
    sv_data_trs.replace('', 0, inplace=True)    
    sv_data_trs.replace('.', 0, inplace=True) 
    sv_data_trs.to_csv("/home/tzhang/trs_test",sep="\t",header=True,index=False)    
    
    # filter low quality sv
    sv_data_trs=sv_data_trs[(sv_data_trs['PR_alt']>=1) & (sv_data_trs['SR_alt']>=1) & ( sv_data_trs['BND_DEPTH'].astype(int) >=1) & ( sv_data_trs['MATE_BND_DEPTH'].astype(int) >=1) & (sv_data_trs['sv_bp_end_POS']!='NAN') & (sv_data_trs['sv_bp_end_POS']!='chrY')]

    # generate sv_info variables (length=bp_st-bp_end; range=ci_st-ci_end; ratio=read_alt/read_ref; diff=read_alt-read_ref)
    sv_data_trs['cipos_range']=sv_data_trs['ci_pos0'].astype(int) -sv_data_trs['ci_pos1'].astype(int)
    sv_data_trs['ciend_range']=sv_data_trs['sv_bp_end_POS'].astype(int) -sv_data_trs['sv_bp_end_END'].astype(int)
     
    sv_data_trs['PR_read_ratio']=abs(sv_data_trs['PR_alt'])/(abs(sv_data_trs['PR_ref'])+abs(sv_data_trs['PR_alt']))
    sv_data_trs['SR_read_ratio']=abs(sv_data_trs['SR_alt'])/(abs(sv_data_trs['SR_ref'])+abs(sv_data_trs['SR_alt']))
    sv_data_trs['read_ratio']=abs(sv_data_trs['PR_alt'])+abs(abs(sv_data_trs['SR_alt']))/(abs(sv_data_trs['SR_ref'])+abs(sv_data_trs['SR_alt'])+abs(sv_data_trs['PR_ref'])+abs(sv_data_trs['PR_alt']))
    sv_data_trs['bnd_depth_ratio']=abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))/(abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))+abs(sv_data_trs['BND_DEPTH'].astype(float)))
    sv_data_trs['PR_read_ratio_diff']=(abs(sv_data_trs['PR_alt'])-abs(sv_data_trs['PR_ref']))/(abs(sv_data_trs['PR_alt'])+abs(sv_data_trs['PR_ref']))
    sv_data_trs['SR_read_ratio_diff']=(abs(sv_data_trs['SR_alt'])-abs(sv_data_trs['SR_ref']))/(abs(sv_data_trs['SR_alt'])+ abs(sv_data_trs['SR_ref']))
    sv_data_trs['read_ratio_diff']=(abs(sv_data_trs['PR_alt'])+ abs(sv_data_trs['SR_alt'])-abs(sv_data_trs['PR_ref'])-abs(sv_data_trs['SR_ref']))/(abs(sv_data_trs['PR_alt'])+ abs(sv_data_trs['SR_alt'])+abs(sv_data_trs['PR_ref'])+abs(sv_data_trs['SR_ref']))
    sv_data_trs['bnd_depth_ratio_diff']=(abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))-abs(sv_data_trs['BND_DEPTH'].astype(float)))/(abs(sv_data_trs['MATE_BND_DEPTH'].astype(float))+abs(sv_data_trs['BND_DEPTH'].astype(float)))
    print(sv_data_trs.shape)   
    
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_trs.shape[0])), k=sv_data_trs.shape[0]) 
    sv_data_trs=sv_data_trs.reset_index(drop=True)
    sv_data_trs=sv_data_trs.reindex(new_index)
    sv_data_trs=sv_data_trs.reset_index(drop=True)     
    # replace na
    sv_data_trs=sv_data_trs.fillna(0)    
    # replace infinite value
    sv_data_trs.replace([np.inf, -np.inf], 0, inplace=True) 
    
    # label benchmark sv 
    sv_data_trs.insert(sv_data_trs.shape[1], "label", list(repeat(0,sv_data_trs.shape[0])), True)
    sv_data_trs.loc[(sv_data_trs['1000_g']=="YVID") | (sv_data_trs['gnomad_g']=="YVID") | (sv_data_trs['dgv_g']=="YVID"),'label']=1   
    sv_data_trs.loc[(sv_data_trs['gnomad_qc']=="YVID") | (sv_data_trs['control_g']=="YVID") | (sv_data_trs['centromere_qc']=="YVID"),'label']=-1  
    sv_data_trs.loc[((sv_data_trs['cytoatlas_s']=="YVID") | (sv_data_trs['cosmic_s']=="YVID")) & (sv_data_trs['control_gall']!="YVID"),'label']=2    
    print(sv_data_trs['label'].value_counts())    
    
    # make summary plot / table
    sv_data_trs['sv_start_bp0']=0
    sv_data_trs['sv_start_bp1']=0
    sv_data_trs['sv_start_bp0']=sv_data_trs['sv_start_bp'].astype(int)-sv_data_trs['ci_pos0'].astype(int)
    sv_data_trs['sv_start_bp1']=sv_data_trs['sv_start_bp'].astype(int)+sv_data_trs['ci_pos1'].astype(int)
    sv_data_trs['sv_chr12']=sv_data_trs['sv_chr1']+":"+sv_data_trs['sv_chr2']
    sv_data_trs.replace('NAN', 0, inplace=True) 
    sv_data_trs.replace('-', 0, inplace=True) 
    sv_data_trs.replace('', 0, inplace=True) 
    
    # drop columns wih rebundant variables
    sv_data_trs_2=sv_data_trs.copy()
    sv_data_trs_2=sv_data_trs_2.loc[:,['sv_type','sv_bp_end_cc1', 'sv_bp_end_cc_v1', 'sv_bp_end_cc_v2',  'sv_bp_end_cc_v3', 'sv_bp_end_cc_v4', 'sv_bp_end_cc_v5', 'sv_bp_end_cc_v6', 'sv_bp_end_cc_v7', 'sv_bp_end_cc_v8', 'sv_bp_end_cc_v9', 'sv_bp_end_cc_v10', 'sv_bp_end_cc_v11',  'sv_bp_end_cc_v12', 'sv_bp_end_cc_v13', 'sv_bp_end_cc_v14',   'sv_bp_end_cc_v15', 'sv_bp_end_cc_v16', 'sv_bp_end_cc_v17', 'sv_bp_end_cc_v18', 'sv_bp_end_cc_v19', 'sv_bp_end_cc_v20', 'sv_bp_end_cc_v21', 'sv_bp_end_cc_v22', 'sv_bp_end_cc_v23', 'sv_bp_end_cc_v24', 'sv_bp_st_cc1', 'sv_bp_st_cc_v1', 'sv_bp_st_cc_v2',  'sv_bp_st_cc_v3', 'sv_bp_st_cc_v4', 'sv_bp_st_cc_v5', 'sv_bp_st_cc_v6',  'sv_bp_st_cc_v7', 'sv_bp_st_cc_v8', 'sv_bp_st_cc_v9', 'sv_bp_st_cc_v10', 'sv_bp_st_cc_v11', 'sv_bp_st_cc_v12', 'sv_bp_st_cc_v13', 'sv_bp_st_cc_v14', 'sv_bp_st_cc_v15', 'sv_bp_st_cc_v16', 'sv_bp_st_cc_v17', 'sv_bp_st_cc_v18', 'sv_bp_st_cc_v19', 'sv_bp_st_cc_v20', 'sv_bp_st_cc_v21', 'sv_bp_st_cc_v22', 'sv_bp_st_cc_v23', 'sv_bp_st_cc_v24', 'cipos_range', 'ciend_range',  'PR_read_ratio', 'SR_read_ratio', 'read_ratio', 'bnd_depth_ratio',  'PR_read_ratio_diff', 'SR_read_ratio_diff', 'read_ratio_diff', 'bnd_depth_ratio_diff', 'label']]
    sv_data_trs_2.columns

#     #  replace false values from data transformation  
#     a=sv_data_trs_2.drop(['sv_type','label'],axis=1).apply(pd.to_numeric, errors="ignore").applymap(lambda x: isinstance(x, float), na_action='ignore')
#     str_cell=[]
#     for i in a.columns:
#         if not a[i].all():
#             for n in range(sv_data_trs_2.shape[0]):
#                 try:
#                     b=float(sv_data_trs_2.loc[n,i])
#                 except:
#                     str_cell.append(str(i)+'&'+str(n))
#     for m in str_cell:
# #        print(m)
#         sv_data_trs_2.loc[int(m.split('&')[1]), m.split('&')[0]]=0   
    
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_trs_2.shape[0])), k=sv_data_trs_2.shape[0]) 
    sv_data_trs_2=sv_data_trs_2.reset_index(drop=True)
    sv_data_trs_2=sv_data_trs_2.reindex(new_index)
    sv_data_trs_2=sv_data_trs_2.reset_index(drop=True)
    return sv_data_trs_2,sv_data_trs

def nontrs_sv_data_transform(sv_data_nontrs, nontrs_sv_cutoff):
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs.shape[0])), k=sv_data_nontrs.shape[0]) 
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)
    sv_data_nontrs=sv_data_nontrs.reindex(new_index)
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)   
    print(sv_data_nontrs.shape) 
    
    # transform database label
    sv_data_nontrs.loc[sv_data_nontrs['1000_gall']=='Not_in_database','1000_gall']=0
    sv_data_nontrs.loc[sv_data_nontrs['1000_g']=='Not_in_database','1000_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['gnomad_gall']=='Not_in_database','gnomad_gall']=0
    sv_data_nontrs.loc[sv_data_nontrs['gnomad_g']=='Not_in_database','gnomad_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['control_gall']=='Not_in_database','control_gall']=0
    sv_data_nontrs.loc[sv_data_nontrs['control_g']=='Not_in_database','control_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['cytoatlas_s']=='Not_in_database','cytoatlas_s']=0
    sv_data_nontrs.loc[sv_data_nontrs['cosmic_s']=='Not_in_database','cosmic_s']=0
    sv_data_nontrs.loc[sv_data_nontrs['gnomad_qc']=='Not_in_database','gnomad_qc']=0
    sv_data_nontrs.loc[sv_data_nontrs['dgv_g']=='Not_in_database','dgv_g']=0
    sv_data_nontrs.loc[sv_data_nontrs['centromere_qc']=='Not_in_database','centromere_qc']=0
    sv_data_nontrs.loc[sv_data_nontrs['uwstl_s']=='Not_in_database','uwstl_s']=0
#     print(sv_data_nontrs['1000_gall'].value_counts())
#     print(sv_data_nontrs['1000_g'].value_counts())
#     print(sv_data_nontrs['gnomad_gall'].value_counts())
#     print(sv_data_nontrs['gnomad_g'].value_counts())
#     print(sv_data_nontrs['control_gall'].value_counts())
#     print(sv_data_nontrs['control_g'].value_counts())
#     print(sv_data_nontrs['cytoatlas_s'].value_counts())
#     print(sv_data_nontrs['cosmic_s'].value_counts())
#     print(sv_data_nontrs['gnomad_qc'].value_counts())
#     print(sv_data_nontrs['centromere_qc'].value_counts())
#     print(sv_data_nontrs['dgv_g'].value_counts())
#     print(sv_data_nontrs['uwstl_s'].value_counts())
#     print(sv_data_nontrs.shape) 
    
    # transform sv_info variables
    sv_data_nontrs['ci_pos0']=0
    sv_data_nontrs['ci_pos1']=0
    sv_data_nontrs['sv_bp_end_POS']=0
    sv_data_nontrs['sv_bp_end_END']=0
    sv_data_nontrs['PR_ref']=0
    sv_data_nontrs['PR_alt']=0
    sv_data_nontrs['SR_ref']=0
    sv_data_nontrs['SR_alt']=0
    sv_data_nontrs['ci_pos0']=sv_data_nontrs['CIPOS'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['ci_pos1']=sv_data_nontrs['CIPOS'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_nontrs['sv_bp_end_POS']=sv_data_nontrs['CIEND'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['sv_bp_end_END']=sv_data_nontrs['CIEND'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
#     print(sv_data_nontrs.shape) 
    
    # read statistics transformation
    sv_data_nontrs.loc[(sv_data_nontrs['RP'].isnull()) | (sv_data_nontrs['RP']=='NA') | (sv_data_nontrs['RP']=='.'),'RP']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['AP'].isnull()) | (sv_data_nontrs['AP']=='NA') | (sv_data_nontrs['AP']=='.'),'AP']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['RS'].isnull()) | (sv_data_nontrs['RS']=='NA') | (sv_data_nontrs['RS']=='.'),'RS']=0          
    sv_data_nontrs.loc[(sv_data_nontrs['AS'].isnull()) | (sv_data_nontrs['AS']=='NA') | (sv_data_nontrs['AS']=='.'),'AS']=0    
    sv_data_nontrs.loc[(sv_data_nontrs['DR'].isnull()) | (sv_data_nontrs['DR']=='NA') | (sv_data_nontrs['DR']=='.'),'DR']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['DV'].isnull()) | (sv_data_nontrs['DV']=='NA') | (sv_data_nontrs['DV']=='.'),'DV']=0  
    sv_data_nontrs.loc[(sv_data_nontrs['RR'].isnull()) | (sv_data_nontrs['RR']=='NA') | (sv_data_nontrs['RR']=='.'),'RR']=0          
    sv_data_nontrs.loc[(sv_data_nontrs['RV'].isnull()) | (sv_data_nontrs['RV']=='NA') | (sv_data_nontrs['RV']=='.'),'RV']=0  
    sv_data_nontrs.to_csv("/home/tzhang/test0",sep="\t",header=True,index=False)     
#     print(sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['AP']!=0),['RP','AP']].shape)
    sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['AP']!=0),'PR']=sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['AP']!=0),['RP','AP']].apply(lambda x: str(x['RP'])+','+str(x['AP']), axis=1)
    sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(',')) | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['AS']!=0),'SR']=sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(',')) | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['AS']!=0),['RS','AS']].apply(lambda x: str(x['RS'])+','+str(x['AS']), axis=1)             
    sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['DV']!=0),'PR']=sv_data_nontrs.loc[((~sv_data_nontrs['PR'].str.contains(',')) | (sv_data_nontrs['PR']=='0,0')) & (sv_data_nontrs['DV']!=0),['DR','DV']].apply(lambda x: str(x['DR'])+','+str(x['DV']), axis=1)
    sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(','))  | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['RV']!=0),'SR']=sv_data_nontrs.loc[((~sv_data_nontrs['SR'].str.contains(','))  | (sv_data_nontrs['SR']=='0,0')) & (sv_data_nontrs['RV']!=0),['RR','RV']].apply(lambda x: str(x['RR'])+','+str(x['RV']), axis=1)  
    sv_data_nontrs.loc[~sv_data_nontrs['PR'].str.contains(','),'PR']='0,0'
    sv_data_nontrs.loc[~sv_data_nontrs['SR'].str.contains(','),'SR']='0,0'       
    sv_data_nontrs.to_csv("/home/tzhang/nontrs_test",sep="\t",header=True,index=False)  
         
    sv_data_nontrs['PR_ref']=sv_data_nontrs['PR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['PR_alt']=sv_data_nontrs['PR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    sv_data_nontrs['SR_ref']=sv_data_nontrs['SR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
    sv_data_nontrs['SR_alt']=sv_data_nontrs['SR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)        

    # replace na
    sv_data_nontrs=sv_data_nontrs.fillna(0)
    sv_data_nontrs.replace('', 0, inplace=True)    
    sv_data_nontrs.replace('.', 0, inplace=True) 
    sv_data_nontrs.to_csv("/home/tzhang/test1",sep="\t",header=True,index=False)    
    
    # filter low quality sv
    sv_data_nontrs['sv_length']=sv_data_nontrs['sv_start_bp'].astype(int)-sv_data_nontrs['sv_end_bp'].astype(int)
    sv_data_nontrs=sv_data_nontrs[sv_data_nontrs['CIEND']!='NA'] 
    #sv_data_nontrs=sv_data_nontrs[(sv_data_nontrs['GT']=='0/1')] 
    #sv_data_nontrs=sv_data_nontrs[sv_data_nontrs['CN']!='NA']
#     print(np.unique(sv_data_nontrs['CN'],return_counts=False))
#     print(np.unique(sv_data_nontrs['PR_alt'],return_counts=False))
#     print(np.unique(sv_data_nontrs['SR_alt'],return_counts=False))
#     print(np.unique(sv_data_nontrs['PR_ref'],return_counts=False))
#     print(np.unique(sv_data_nontrs['SR_ref'],return_counts=False))
#     print(np.unique(sv_data_nontrs['SUPP'],return_counts=False))      
#     print(sv_data_nontrs['CN'].value_counts())
#     print(sv_data_nontrs['PR_alt'].value_counts())
#     print(sv_data_nontrs['SR_alt'].value_counts())
#     print(sv_data_nontrs['PR_ref'].value_counts())
#     print(sv_data_nontrs['SR_ref'].value_counts())
#     print(sv_data_nontrs['SUPP'].value_counts())    
    sv_data_nontrs=sv_data_nontrs[ (sv_data_nontrs['CN']!='2') | (abs(sv_data_nontrs['PR_alt'].astype(int))>1) | (abs(sv_data_nontrs['SR_alt'].astype(int))>0) | (abs(sv_data_nontrs['PR_ref'].astype(int))>0) |  (abs(sv_data_nontrs['SR_ref'].astype(int))>0) | (sv_data_nontrs['SUPP'].astype(int)>1)]
   
    # generate sv_info variables (length=bp_st-bp_end; range=ci_st-ci_end; ratio=read_alt/read_ref; diff=read_alt-read_ref)
    sv_data_nontrs['cipos_range']=sv_data_nontrs['ci_pos0'].astype(int) -sv_data_nontrs['ci_pos1'].astype(int)
    sv_data_nontrs['ciend_range']=sv_data_nontrs['sv_bp_end_POS'].astype(int) -sv_data_nontrs['sv_bp_end_END'].astype(int)
     
    sv_data_nontrs['PR_read_ratio']=abs(sv_data_nontrs['PR_alt'])/(abs(sv_data_nontrs['PR_ref'])+abs(sv_data_nontrs['PR_alt']))
    sv_data_nontrs['SR_read_ratio']=abs(sv_data_nontrs['SR_alt'])/(abs(sv_data_nontrs['SR_ref'])+abs(sv_data_nontrs['SR_alt']))
    sv_data_nontrs['read_ratio']=abs(sv_data_nontrs['PR_alt'])+abs(abs(sv_data_nontrs['SR_alt']))/(abs(sv_data_nontrs['SR_ref'])+abs(sv_data_nontrs['SR_alt'])+abs(sv_data_nontrs['PR_ref'])+abs(sv_data_nontrs['PR_alt']))
    sv_data_nontrs['PR_read_ratio_diff']=(abs(sv_data_nontrs['PR_alt'])-abs(sv_data_nontrs['PR_ref']))/(abs(sv_data_nontrs['PR_alt'])+abs(sv_data_nontrs['PR_ref']))
    sv_data_nontrs['SR_read_ratio_diff']=(abs(sv_data_nontrs['SR_alt'])-abs(sv_data_nontrs['SR_ref']))/(abs(sv_data_nontrs['SR_alt'])+ abs(sv_data_nontrs['SR_ref']))
    sv_data_nontrs['read_ratio_diff']=(abs(sv_data_nontrs['PR_alt'])+ abs(sv_data_nontrs['SR_alt'])-abs(sv_data_nontrs['PR_ref'])-abs(sv_data_nontrs['SR_ref']))/(abs(sv_data_nontrs['PR_alt'])+ abs(sv_data_nontrs['SR_alt'])+abs(sv_data_nontrs['PR_ref'])+abs(sv_data_nontrs['SR_ref']))
    
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs.shape[0])), k=sv_data_nontrs.shape[0]) 
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)
    sv_data_nontrs=sv_data_nontrs.reindex(new_index)
    sv_data_nontrs=sv_data_nontrs.reset_index(drop=True)       
    # replace na
    sv_data_nontrs=sv_data_nontrs.fillna(0)
    # replace infinite value
    sv_data_nontrs.replace([np.inf, -np.inf], 0, inplace=True)
    sv_data_nontrs.to_csv("/home/tzhang/test2"+str(nontrs_sv_cutoff)+"bf",sep="\t",header=True,index=False)   
    print(sv_data_nontrs.shape)  
        
    # label benchmark sv 
    print(sv_data_nontrs.shape)      
    sv_data_nontrs=sv_data_nontrs[abs(sv_data_nontrs['sv_length'])>=100000]    
    print(sv_data_nontrs.shape)      
    sv_data_nontrs=sv_data_nontrs[(sv_data_nontrs['GT']=='0/1')]
    print(sv_data_nontrs.shape)     
    sv_data_nontrs.insert(sv_data_nontrs.shape[1], "label", list(repeat(0,sv_data_nontrs.shape[0])), True)
    sv_data_nontrs.loc[(sv_data_nontrs['1000_g'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['gnomad_g'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['dgv_g'].astype(float)>nontrs_sv_cutoff),'label']=1
    sv_data_nontrs.loc[( (sv_data_nontrs['cosmic_s'].astype(float)>nontrs_sv_cutoff))  & (sv_data_nontrs['dgv_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['1000_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['gnomad_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['control_g'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['1000_gall'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['gnomad_gall'].astype(float)<nontrs_sv_cutoff) & (sv_data_nontrs['control_gall'].astype(float)<nontrs_sv_cutoff) ,'label']=2    
    sv_data_nontrs.loc[(sv_data_nontrs['gnomad_qc'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['control_g'].astype(float)>nontrs_sv_cutoff) | (sv_data_nontrs['centromere_qc'].astype(float)>nontrs_sv_cutoff),'label']=-1
    sv_data_nontrs.to_csv("/home/tzhang/test2"+str(nontrs_sv_cutoff)+"af",sep="\t",header=True,index=False)
    
    print(sv_data_nontrs['label'].value_counts())      
    sv_data_nontrs=sv_data_nontrs.rename(columns={'QUAL':'svtyper_score','SUPP':'sv_caller_supp'})  
    sv_data_nontrs.replace('NAN', 0, inplace=True) 
    sv_data_nontrs.replace('-', 0, inplace=True)  
    sv_data_nontrs.replace('', 0, inplace=True)  
    
    # drop columns wih rebundant variables
    sv_data_nontrs_2=sv_data_nontrs.copy()
    sv_data_nontrs_2=sv_data_nontrs_2.loc[:,['sv_type', 'svtyper_score', 'sv_bp_end_cc1', 'sv_bp_end_cc_v1','sv_bp_end_cc_v2', 'sv_bp_end_cc_v3', 'sv_bp_end_cc_v4','sv_bp_end_cc_v5', 'sv_bp_end_cc_v6', 'sv_bp_end_cc_v7', 'sv_bp_end_cc_v8', 'sv_bp_end_cc_v9', 'sv_bp_end_cc_v10', 'sv_bp_end_cc_v11', 'sv_bp_end_cc_v12', 'sv_bp_end_cc_v13', 'sv_bp_end_cc_v14', 'sv_bp_end_cc_v15', 'sv_bp_end_cc_v16','sv_bp_end_cc_v17', 'sv_bp_end_cc_v18', 'sv_bp_end_cc_v19', 'sv_bp_end_cc_v20', 'sv_bp_end_cc_v21', 'sv_bp_end_cc_v22', 'sv_bp_end_cc_v23', 'sv_bp_end_cc_v24', 'sv_bp_st_cc1', 'sv_bp_st_cc_v1', 'sv_bp_st_cc_v2', 'sv_bp_st_cc_v3', 'sv_bp_st_cc_v4', 'sv_bp_st_cc_v5', 'sv_bp_st_cc_v6', 'sv_bp_st_cc_v7', 'sv_bp_st_cc_v8', 'sv_bp_st_cc_v9', 'sv_bp_st_cc_v10', 'sv_bp_st_cc_v11', 'sv_bp_st_cc_v12', 'sv_bp_st_cc_v13', 'sv_bp_st_cc_v14', 'sv_bp_st_cc_v15', 'sv_bp_st_cc_v16', 'sv_bp_st_cc_v17', 'sv_bp_st_cc_v18', 'sv_bp_st_cc_v19', 'sv_bp_st_cc_v20', 'sv_bp_st_cc_v21', 'sv_bp_st_cc_v22', 'sv_bp_st_cc_v23', 'sv_bp_st_cc_v24', 'sv_caller_supp', 'cipos_range', 'ciend_range', 'PR_read_ratio', 'SR_read_ratio', 'read_ratio', 'PR_read_ratio_diff',  'SR_read_ratio_diff', 'read_ratio_diff', 'label']]
    sv_data_nontrs_2.columns
    
    # re-arrange the data 
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs_2.shape[0])), k=sv_data_nontrs_2.shape[0]) 
    sv_data_nontrs_2=sv_data_nontrs_2.reset_index(drop=True)
    sv_data_nontrs_2=sv_data_nontrs_2.reindex(new_index)
    sv_data_nontrs_2=sv_data_nontrs_2.reset_index(drop=True)    
    #sv_data_nontrs_2=sv_data_nontrs_2[sv_data_nontrs_2['sv_type']!='INV']
    #sv_data_nontrs_2.loc[sv_data_nontrs_2['svtyper_score']=='.','svtyper_score']=0
    sv_data_nontrs_2.to_csv("/home/tzhang/test4",sep="\t",header=True,index=False)       
    return sv_data_nontrs_2,sv_data_nontrs

def sv_data_summary_plot(sv_data):
    # count labels by  sv type
    sv_datat2=sv_data.copy()
    sv_datat2.insert(sv_datat2.shape[1], "count", list(repeat(1,sv_datat2.shape[0])), True)
    sv_datat2=sv_datat2.loc[:,['sv_type','label','count']]
    print(sv_datat2.groupby(['sv_type','label']).sum())
    SA=sv_datat2.loc[sv_datat2['label']==-1,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    UL=sv_datat2.loc[sv_datat2['label']==0,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    TG=sv_datat2.loc[sv_datat2['label']==1,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    TC=sv_datat2.loc[sv_datat2['label']==2,['sv_type','count']].groupby(['sv_type']).sum()['count'].values
    index=sv_datat2['sv_type'].unique()
    df1 = pd.DataFrame({'class -1': SA,'class 1': TG,'class 2': TC}, index=index)
    print(df1.iloc[0,:].sum())    
    ax = df1.plot.bar(stacked=True, figsize=(5, 7),width = 0.2,legend=None)
    print(df1)
    for idx, sv in enumerate(df1.index):
        df2= df1[df1.index==sv]
        d0=0
        for i in range(df2.shape[1]):
            data=df2.iloc[0,i]
            ax.text(x=idx-0.1 , y = d0+data/2*pow(-2,i), s=f"{data}" , fontdict=dict(fontsize=10))
            d0=d0+data
    ax.set_yscale('linear')
    ax.ticklabel_format(style='plain', axis='y')
    ax.set_ylabel("SV number")
    ax.set_ylim([0,pow(10,len(str(df1.iloc[0,:].sum())))*1.1])
    fig = ax.figure 
    return fig
