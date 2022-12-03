import os
import sys
import pandas as pd
import numpy as np
import pathlib
import snakemake.io
from snakemake.utils import validate
from typing import Dict, Union, List

configfile: "config.yaml"
    
cohort_name = config['cohort_name']  
samples_information = pd.read_csv(config['sample_list'], sep='\t', header=None,index_col=False)
samples_information.columns=['id','sex']
SAMPLES = list(samples_information['id'])
GENDERS = list(samples_information['sex'])
SAMPLES_vector='@'.join(str(sm) for sm in SAMPLES)

MAIN_DIR = config['main_dir']
INPUT_DIR = config['main_dir']+'/in'
OUTPUT_DIR = config['main_dir']+'/out'
LOG_DIR = config['main_dir']+'/out/log'
CYTO_SV_ML_DIR = config['cyto_sv_ml_dir']
SOFTWARE_DIR = config['cyto_sv_ml_dir']+'/software'
DATABASE_DIR = config['cyto_sv_ml_dir']+'/SV_database'
parliment_docker = config['parliment_docker']
chromoseq_docker = config['chromoseq_docker']
parliment2_sv_callers = config['parliment2_sv_callers']
chromoseq_sv_callers = config['chromoseq_sv_callers']
all_callers=chromoseq_sv_callers+parliment2_sv_callers
all_callers_svtyper=['manta', 'delly', 'cnvnator', 'breakdancer']
size=int(config['size'])
SIZE_K=round(size/1000)
#report: OUTPUT_DIR+"/report/workflow.rst"

rule all:
    input:
        expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_{sv_type}_sv_ml_metrics_sub.csv", cohort_name=cohort_name, sv_type=['trs','nontrs'])
        
#  checkpoint for all sample sv data   
checkpoint all_sample_sv_ready:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES, size_k=SIZE_K)  
    output: 
        pathlib.Path(OUTPUT_DIR+"/log_files/sample_sv_ready.out")
    run: 
        shell('echo {SAMPLES_vector} >> {output}')
        
def check_sample_file(*wildcards):
     return checkpoints.all_sample_sv_ready.get().output        

# combine all sample sv           
rule all_sample_sv_combine:
    input:
        check_sample_file
    output:
        expand(OUTPUT_DIR+"/{cohort_name}/{cohort_name}.sv.all.combine_all", cohort_name=cohort_name)        
    shell:
        """  
        cat {input} && bash {CYTO_SV_ML_DIR}/Pipeline_script/all_sample_sv_combine.sh {MAIN_DIR} {cohort_name} {input} 
        """               
               
# run cyto-sv-ml model     
rule cyto_sv_ml:
    input:
        expand(OUTPUT_DIR+"/{cohort_name}/{cohort_name}.sv.all.combine_all", cohort_name=cohort_name)   
    output:
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_sv_summary_plot.pdf", cohort_name=cohort_name,k=0), category="sv data summary", subcategory="data",labels={"data name" : "sv type distribution","sv type": "trs","data type": "plot" }),    
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_sv_summary_plot.pdf", cohort_name=cohort_name,k=0), category="sv data summary", subcategory="data",labels={"data name" : "sv type distribution","sv type": "nontrs","data type": "plot" }),         
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_sv_ml_metrics_sub.csv", cohort_name=cohort_name), category="sv model summary", subcategory="model", labels={"data name" : "model performance metrics","sv type":'trs', "data type": "table" }), 
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_sv_ml_metrics_sub.csv", cohort_name=cohort_name), category="sv model summary", subcategory="model", labels={"data name" : "model performance metrics","sv type":'nontrs', "data type": "table" }), 
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_{k}_ts_EXP/learner_fold_0_shap_summary.png", cohort_name=cohort_name,k=0), category="sv model summary", subcategory="model",labels={"data name" : "shap feature importance","sv type": "trs","data type": "plot" }),    
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_{k}_ts_EXP/learner_fold_0_shap_summary.png", cohort_name=cohort_name,k=0), category="sv model summary", subcategory="model",labels={"data name" : "shap feature importance","sv type": "nontrs","data type": "plot" }),        
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_{k}_ts_model_confusion_matrix.pdf", cohort_name=cohort_name,k=0), category="sv model summary", subcategory="model",labels={"data name" : "model confusion matrix ","sv type": "trs","data type": "plot" }),    
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_{k}_ts_model_confusion_matrix.pdf", cohort_name=cohort_name,k=0), category="sv model summary", subcategory="model",labels={"data name" : "model confusion matrix","sv type": "nontrs","data type": "plot" }),
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_trs_{k}_ts_model_aucroc_curve.pdf", cohort_name=cohort_name,k=0), category="sv model summary", subcategory="model",labels={"data name" : "model_aucroc_curve","sv type": "trs","data type": "plot" }),    
        report(expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_nontrs_{k}_ts_model_aucroc_curve.pdf", cohort_name=cohort_name,k=0), category="sv model summary", subcategory="model",labels={"data name" : "model aucroc curve","sv type": "nontrs","data type": "plot" })         
    params:
        kfs=config['kfolds_shuffle'] 
    shell:
        """
        sudo mkdir -p {OUTPUT_DIR}/{cohort_name}/cyto_sv_ml &&
        python {CYTO_SV_ML_DIR}/Pipeline_script/CYTO-SV-Auto-ML.py -s {cohort_name} -o {OUTPUT_DIR}/{cohort_name} -k {params.kfs} 
        """             
