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
    run: 
        shell('echo {input} >> {output}')

def check_sample_file(*wildcards):
     return checkpoints.all_sample_sv_ready.get().output        

# combine all sample sv           
rule all_sv_combine:
    input:
        check_sample_file, 
        sv_all_combine=expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES) 
    output:
        expand(OUTPUT_DIR+"/{sample_all}.sv.all.combine_all", sample_all=sample_all)       
    shell:
        """        
         cat {input.sv_all_combine} | awk '(FNR==1)||($1!~"sv_id"){print $0}' >> {output}
        """               
               
# run cyto-sv-ml model     
rule cyto-sv-ml:
#    conda:
#        "conda-py37.yaml"
    output:
       report(expand(OUTPUT_DIR+"/{sample}/{sample}.model_{analyses}.pdf", sample=SAMPLES, analyses=["type_class_summary","score_metrics", "confusion_matrix", "aucroc_curve"]))
    shell:
        'python {main_dir}/software/CYTO-SV-ML/Pipeline_script/CYTO-SV-Auto-ML_tuning.py -s {sample_all} -o {sample_all}_ts -n 10'           
    
