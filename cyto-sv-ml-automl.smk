import os
import sys
import pandas as pd
import numpy as np
import pathlib
import snakemake.io
from snakemake.utils import validate
from typing import Dict, Union, List

configfile: "config.yaml"
    
sample_all = config['cohort_name']  
samples_information = pd.read_csv(config['sample_list'], sep='\t', header=None,index_col=False)
samples_information.columns=['id','sex']
SAMPLES = list(samples_information['id'])
GENDERS = list(samples_information['sex'])
SAMPLES_vector='@'.join(str(sm) for sm in SAMPLES)

#  checkpoint for all sample sv data   
checkpoint all_sample_sv_ready:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES)  
    output: 
        pathlib.Path(OUTPUT_DIR+"/log_files/sample_sv_ready.out")
    para:
        SAMPLES_vector=SAMPLES_vector
    run: 
        shell('echo {input} >> {output}')

def check_sample_file(*wildcards):
     return checkpoints.all_sample_sv_ready.get().output        

# combine all sample sv           
rule all_sample_sv_combine:
    input:
        check_sample_file, 
        sv_all_combine=expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES) 
    output:
        expand(OUTPUT_DIR+"/{sample_all}.sv.all.combine_all", sample_all=sample_all)       
    shell:
        """        
         bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/all_sample_sv_combine.sh {main_dir} {sample_all}
        """               
               
# run cyto-sv-ml model     
rule cyto-sv-ml:
#    conda:
#        "conda-py37.yaml"
    output:
       report(expand(OUTPUT_DIR+"/"+sample_all+"_{sv_type}_sv_ml_metrics_sub.csv", sv_type=['TRS','NONTRS']),
              expand(OUTPUT_DIR+"/"+sample_all+"_{sv_type}_svtype_class_summary.pdf", sv_type=['TRS','NONTRS']),
              expand(OUTPUT_DIR+"/"+sample_all+"_{sv_type}_{analyses}.pdf", sv_type=['TRS','NONTRS'], analyses=["confusion_matrix", "aucroc_curve"]))
    shell:
        'python {main_dir}/software/CYTO-SV-ML/Pipeline_script/CYTO-SV-Auto-ML_tuning.py -s {sample_all} -o {sample_all}_ts -n 10'           
    
