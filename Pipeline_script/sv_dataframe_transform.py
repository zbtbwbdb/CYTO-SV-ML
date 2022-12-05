import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random 
from itertools import repeat

def trs_sv_data_transform(sv_data_trs):
    # transform database label
    sv_db_list=['1000_gall','1000_g','gnomad_gall','gnomad_g','control_gall', 'control_g', 'cytoatlas_s','cosmic_s','gnomad_qc', 'dgv_g', 'centromere_qc','uwstl_s' ]
    for sv_db in sv_db_list:
        sv_data_trs.loc[sv_data_trs[sv_db]=='NAN',sv_db]="Not_in_database"
    sv_data_trs.loc[sv_data_trs['1000_gall']!='Not_in_database','1000_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['1000_g']!='Not_in_database','1000_g']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_gall']!='Not_in_database','gnomad_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_g']!='Not_in_database','gnomad_g']="YVID"
    sv_data_trs.loc[sv_data_trs['control_gall']!='Not_in_database','control_gall']="YVID"
    sv_data_trs.loc[sv_data_trs['control_g']!='Not_in_database','control_g']="YVID"
    sv_data_trs.loc[sv_data_trs['cytoatlas_s']!='Not_in_database','cytoatlas_s']="YVID"
    sv_data_trs.loc[sv_data_trs['cosmic_s']!='Not_in_database','cosmic_s']="YVID"
    sv_data_trs.loc[sv_data_trs['gnomad_qc']!='Not_in_database','gnomad_qc']="YVID"
    sv_data_trs.loc[sv_data_trs['dgv_g']!='Not_in_database','dgv_g']="YVID"
    sv_data_trs.loc[sv_data_trs['centromere_qc']!='Not_in_database','centromere_qc']="YVID"
    sv_data_trs.loc[sv_data_trs['uwstl_s']!='Not_in_database','uwstl_s']="YVID"
    
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
    if 'RP' in sv_data_trs.columns:
        sv_data_trs.loc[sv_data_trs['AS']=='NA' ,'AS']=0         
        sv_data_trs['PR_ref']=sv_data_trs['RP'].astype(int)
        sv_data_trs['PR_alt']=sv_data_trs['AP'].astype(int)
        sv_data_trs['SR_ref']=sv_data_trs['RS'].astype(int)
        sv_data_trs['SR_alt']=sv_data_trs['AS'].astype(int) 
    else:
        sv_data_trs.loc[sv_data_trs['SR']=='NA' ,'SR']='0,0'          
        sv_data_trs['PR_ref']=sv_data_trs['PR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
        sv_data_trs['PR_alt']=sv_data_trs['PR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
        sv_data_trs['SR_ref']=sv_data_trs['SR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
        sv_data_trs['SR_alt']=sv_data_trs['SR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    
    # filter low quality sv
    sv_data_trs=sv_data_trs[(sv_data_trs['PR_alt']>1) & (sv_data_trs['SR_alt']>1) & ( sv_data_trs['BND_DEPTH'].astype(int) > 1) & ( sv_data_trs['MATE_BND_DEPTH'].astype(int) > 1) & (sv_data_trs['sv_bp_end_POS']!='NAN') & (sv_data_trs['sv_bp_end_POS']!='chrY')]
    
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
    
    # replace na
    sv_data_trs=sv_data_trs.fillna(0)
    # replace infinite value
    sv_data_trs.replace([np.inf, -np.inf], 0, inplace=True)
    
    # label benchmark sv 
    sv_data_trs.insert(sv_data_trs.shape[1], "label", list(repeat(0,sv_data_trs.shape[0])), True)
    sv_data_trs.loc[(sv_data_trs['gnomad_qc']=="YVID") | (sv_data_trs['control_g']=="YVID") | (sv_data_trs['centromere_qc']=="YVID"),'label']=-1
    sv_data_trs.loc[(sv_data_trs['1000_g']=="YVID") | (sv_data_trs['gnomad_g']=="YVID") | (sv_data_trs['dgv_g']=="YVID") | (sv_data_trs['1000_gall']=="YVID") | (sv_data_trs['gnomad_gall']=="YVID"),'label']=1
    sv_data_trs.loc[((sv_data_trs['cytoatlas_s']=="YVID") | (sv_data_trs['cosmic_s']=="YVID")) & (sv_data_trs['1000_gall']!="YVID") & (sv_data_trs['gnomad_gall']!="YVID") & (sv_data_trs['control_gall']!="YVID") ,'label']=2
    
    
    # make summary plot / table
    sv_data_trs['sv_start_bp0']=0
    sv_data_trs['sv_start_bp1']=0
    sv_data_trs['sv_start_bp0']=sv_data_trs['sv_start_bp'].astype(int)-sv_data_trs['ci_pos0'].astype(int)
    sv_data_trs['sv_start_bp1']=sv_data_trs['sv_start_bp'].astype(int)+sv_data_trs['ci_pos1'].astype(int)
    sv_data_trs['sv_chr12']=sv_data_trs['sv_chr1']+":"+sv_data_trs['sv_chr2']
    
    # drop columns wih rebundant variables
    sv_data_trs_2=sv_data_trs.copy()
    sv_data_trs_2=sv_data_trs_2.loc[:,['sv_type','sv_bp_end_cc1', 'sv_bp_end_cc_v1', 'sv_bp_end_cc_v2',  'sv_bp_end_cc_v3', 'sv_bp_end_cc_v4', 'sv_bp_end_cc_v5', 'sv_bp_end_cc_v6', 'sv_bp_end_cc_v7', 'sv_bp_end_cc_v8', 'sv_bp_end_cc_v9', 'sv_bp_end_cc_v10', 'sv_bp_end_cc_v11',  'sv_bp_end_cc_v12', 'sv_bp_end_cc_v13', 'sv_bp_end_cc_v14',   'sv_bp_end_cc_v15', 'sv_bp_end_cc_v16', 'sv_bp_end_cc_v17', 'sv_bp_end_cc_v18', 'sv_bp_end_cc_v19', 'sv_bp_end_cc_v20', 'sv_bp_end_cc_v21', 'sv_bp_end_cc_v22', 'sv_bp_end_cc_v23', 'sv_bp_end_cc_v24', 'sv_bp_st_cc1', 'sv_bp_st_cc_v1', 'sv_bp_st_cc_v2',  'sv_bp_st_cc_v3', 'sv_bp_st_cc_v4', 'sv_bp_st_cc_v5', 'sv_bp_st_cc_v6',  'sv_bp_st_cc_v7', 'sv_bp_st_cc_v8', 'sv_bp_st_cc_v9', 'sv_bp_st_cc_v10', 'sv_bp_st_cc_v11', 'sv_bp_st_cc_v12', 'sv_bp_st_cc_v13', 'sv_bp_st_cc_v14', 'sv_bp_st_cc_v15', 'sv_bp_st_cc_v16', 'sv_bp_st_cc_v17', 'sv_bp_st_cc_v18', 'sv_bp_st_cc_v19', 'sv_bp_st_cc_v20', 'sv_bp_st_cc_v21', 'sv_bp_st_cc_v22', 'sv_bp_st_cc_v23', 'sv_bp_st_cc_v24', 'cipos_range', 'ciend_range',  'PR_read_ratio', 'SR_read_ratio', 'read_ratio', 'bnd_depth_ratio',  'PR_read_ratio_diff', 'SR_read_ratio_diff', 'read_ratio_diff', 'bnd_depth_ratio_diff', 'label']]
    sv_data_trs_2.columns
    
    # re-arrange the data
    random.seed(0)
    new_index=random.sample(list(range(sv_data_trs_2.shape[0])), k=sv_data_trs_2.shape[0]) 
    sv_data_trs_2=sv_data_trs_2.reset_index(drop=True)
    sv_data_trs_2=sv_data_trs_2.reindex(new_index)
    sv_data_trs_2=sv_data_trs_2.reset_index(drop=True)
    return sv_data_trs_2

def nontrs_sv_data_transform(sv_data_nontrs):
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
    if 'RP' in sv_data_nontrs.columns:    
        sv_data_nontrs.loc[sv_data_nontrs['AS']=='NA' ,'AS']=0        
        sv_data_nontrs['PR_ref']=sv_data_nontrs['RP'].astype(int)
        sv_data_nontrs['PR_alt']=sv_data_nontrs['AP'].astype(int)
        sv_data_nontrs['SR_ref']=sv_data_nontrs['RS'].astype(int)
        sv_data_nontrs['SR_alt']=sv_data_nontrs['AS'].astype(int)  
    else:
        sv_data_nontrs.loc[sv_data_nontrs['SR']=='NA' ,'SR']='0,0'          
        sv_data_nontrs['PR_ref']=sv_data_nontrs['PR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
        sv_data_nontrs['PR_alt']=sv_data_nontrs['PR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
        sv_data_nontrs['SR_ref']=sv_data_nontrs['SR'].apply(lambda x: int(x.split(',')[0]) if ',' in x else 0)
        sv_data_nontrs['SR_alt']=sv_data_nontrs['SR'].apply(lambda x: int(x.split(',')[1]) if ',' in x else 0)
    
    # filter low quality sv
    sv_data_nontrs['sv_length']=sv_data_nontrs['sv_start_bp'].astype(int)-sv_data_nontrs['sv_end_bp'].astype(int)
    sv_data_nontrs=sv_data_nontrs[(abs(sv_data_nontrs['sv_length'])>=100000) ] 
    sv_data_nontrs=sv_data_nontrs[sv_data_nontrs['CIEND']!='NA'] 
    #sv_data_nontrs=sv_data_nontrs[(sv_data_nontrs['GT']=='0/1')] 
    #sv_data_nontrs=sv_data_nontrs[sv_data_nontrs['CN']!='NA']
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
    
    # replace na
    sv_data_nontrs=sv_data_nontrs.fillna(0)
    # replace infinite value
    sv_data_nontrs.replace([np.inf, -np.inf], 0, inplace=True)
    
    # label benchmark sv 
    sv_data_nontrs.insert(sv_data_nontrs.shape[1], "label", list(repeat(0,sv_data_nontrs.shape[0])), True)
    sv_data_nontrs.loc[(sv_data_nontrs['gnomad_qc'].astype(float)>0.9) | (sv_data_nontrs['control_g'].astype(float)>0.9) | (sv_data_nontrs['centromere_qc'].astype(float)>0.9),'label']=-1
    sv_data_nontrs.loc[(sv_data_nontrs['1000_g'].astype(float)>0.7) | (sv_data_nontrs['gnomad_g'].astype(float)>0.7) | (sv_data_nontrs['dgv_g'].astype(float)>0.7) | (sv_data_nontrs['1000_gall'].astype(float)>0.7) | (sv_data_nontrs['gnomad_gall'].astype(float)>0.7),'label']=1
    sv_data_nontrs.loc[( (sv_data_nontrs['cytoatlas_s'].astype(float)>0.9) | (sv_data_nontrs['cosmic_s'].astype(float)>0.9)) & (sv_data_nontrs['1000_gall'].astype(float)<0.9) & (sv_data_nontrs['gnomad_gall'].astype(float)<0.9) & (sv_data_nontrs['control_gall'].astype(float)<0.9) ,'label']=2
    
    # drop columns wih rebundant variables
    sv_data_nontrs_2=sv_data_nontrs.copy()
    sv_data_nontrs_2=sv_data_nontrs_2.rename(columns={'QUAL.1':'svtyper_score','SUPP':'sv_caller_supp'})    
    sv_data_nontrs_2=sv_data_nontrs_2.loc[:,['sv_type', 'svtyper_score', 'sv_bp_end_cc1', 'sv_bp_end_cc_v1','sv_bp_end_cc_v2', 'sv_bp_end_cc_v3', 'sv_bp_end_cc_v4','sv_bp_end_cc_v5', 'sv_bp_end_cc_v6', 'sv_bp_end_cc_v7', 'sv_bp_end_cc_v8', 'sv_bp_end_cc_v9', 'sv_bp_end_cc_v10', 'sv_bp_end_cc_v11', 'sv_bp_end_cc_v12', 'sv_bp_end_cc_v13', 'sv_bp_end_cc_v14', 'sv_bp_end_cc_v15', 'sv_bp_end_cc_v16','sv_bp_end_cc_v17', 'sv_bp_end_cc_v18', 'sv_bp_end_cc_v19', 'sv_bp_end_cc_v20', 'sv_bp_end_cc_v21', 'sv_bp_end_cc_v22', 'sv_bp_end_cc_v23', 'sv_bp_end_cc_v24', 'sv_bp_st_cc1', 'sv_bp_st_cc_v1', 'sv_bp_st_cc_v2', 'sv_bp_st_cc_v3', 'sv_bp_st_cc_v4', 'sv_bp_st_cc_v5', 'sv_bp_st_cc_v6', 'sv_bp_st_cc_v7', 'sv_bp_st_cc_v8', 'sv_bp_st_cc_v9', 'sv_bp_st_cc_v10', 'sv_bp_st_cc_v11', 'sv_bp_st_cc_v12', 'sv_bp_st_cc_v13', 'sv_bp_st_cc_v14', 'sv_bp_st_cc_v15', 'sv_bp_st_cc_v16', 'sv_bp_st_cc_v17', 'sv_bp_st_cc_v18', 'sv_bp_st_cc_v19', 'sv_bp_st_cc_v20', 'sv_bp_st_cc_v21', 'sv_bp_st_cc_v22', 'sv_bp_st_cc_v23', 'sv_bp_st_cc_v24', 'sv_caller_supp', 'cipos_range', 'ciend_range', 'PR_read_ratio', 'SR_read_ratio', 'read_ratio', 'PR_read_ratio_diff',  'SR_read_ratio_diff', 'read_ratio_diff', 'label']]
    sv_data_nontrs_2.columns
    
    # re-arrange the data 
    random.seed(0)
    new_index=random.sample(list(range(sv_data_nontrs_2.shape[0])), k=sv_data_nontrs_2.shape[0]) 
    sv_data_nontrs_2=sv_data_nontrs_2.reset_index(drop=True)
    sv_data_nontrs_2=sv_data_nontrs_2.reindex(new_index)
    sv_data_nontrs_2=sv_data_nontrs_2.reset_index(drop=True)    
    sv_data_nontrs_2=sv_data_nontrs_2[sv_data_nontrs_2['sv_type']!='INV']
    sv_data_nontrs_2.loc[sv_data_nontrs_2['svtyper_score']=='.','svtyper_score']=0
    return sv_data_nontrs_2

def sv_data_summary_plot(sv_data):
    # count labels by  sv type
    sv_datat2=sv_data.copy()
    sv_datat2.insert(sv_datat2.shape[1], "count", list(repeat(1,sv_datat2.shape[0])), True)
    sv_datat2=sv_datat2.loc[:,['sv_type','label','count']]
    print(sv_datat2.groupby(['label']).sum())
    SA=sv_datat2.loc[sv_datat2['label']==-1,['count']].sum().values
    UL=sv_datat2.loc[sv_datat2['label']==0,['count']].sum().values
    TG=sv_datat2.loc[sv_datat2['label']==1,['count']].sum().values
    TC=sv_datat2.loc[sv_datat2['label']==2,['count']].sum().values
    index=sv_datat2['sv_type'].unique()
    df1 = pd.DataFrame({'class -1': SA,'class 1': TG,'class 2': TC}, index=index)
    print(df1.iloc[0,:].sum())    
    ax = df1.plot.bar(stacked=True, figsize=(5, 7),width = 0.2,legend=None)
    d0=0
    for i in range(df1.shape[1]):
        data=df1.iloc[0,i]
        ax.text(x=-0.05, y = d0+data/2 , s=f"{data}" , fontdict=dict(fontsize=10))
        d0=d0+data
    ax.set_yscale('linear')
    ax.ticklabel_format(style='plain', axis='y')
    ax.set_ylabel("SV number")
    ax.set_ylim([0,pow(10,len(str(df1.iloc[0,:].sum())))])
    fig = ax.figure
    return fig
