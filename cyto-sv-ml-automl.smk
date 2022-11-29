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

#INPUT_DIR = pathlib.Path(config['main_dir']+'/in')
#OUTPUT_DIR = pathlib.Path(config['main_dir']+'/out')
main_dir = config['main_dir']
INPUT_DIR = config['main_dir']+'/in'
OUTPUT_DIR = config['main_dir']+'/out'
LOG_DIR = config['main_dir']+'/out/log'
REF_DIR = config['main_dir']+'/reference'
SOFTWARE_DIR = config['main_dir']+'/software'
DATABASE_DIR = config['main_dir']+'/SV_database'
parliment2_sv_callers = config['parliment2_sv_callers']
chromoseq_sv_callers = config['chromoseq_sv_callers']
all_callers=chromoseq_sv_callers+parliment2_sv_callers
all_callers_svtyper=['manta', 'delly', 'cnvnator', 'breakdancer']
# print(os.path.join(OUTPUT_DIR,"/log_files/sample_sv_ready.out"))
# print(pathlib.Path(OUTPUT_DIR+"/log_files/sample_sv_ready.out"))
size=int(config['size'])
#report: OUTPUT_DIR+"/report/workflow.rst"

rule all:
    input:
        expand(OUTPUT_DIR+"/"+sample_all+"_{sv_type}_sv_ml_metrics_sub.csv", sv_type=['TRS','NONTRS'])  
        
#  checkpoint for all sample sv data   
checkpoint all_sample_sv_ready:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES)  
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
        expand(OUTPUT_DIR+"/{sample_all}.sv.all.combine_all", sample_all=sample_all)        
    shell:
        """  
        cat {input} && bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/all_sample_sv_combine.sh {main_dir} {sample_all} {input} 
        """               
               
# run cyto-sv-ml model     
rule cyto_sv_ml:
    input:
        expand(OUTPUT_DIR+"/{sample_all}.sv.all.combine_all", sample_all=sample_all)   
    params:
        py39_dir=config['py39_dir']
    output:
        report(expand(OUTPUT_DIR+"/"+sample_all+"_{sv_type}_sv_ml_metrics_sub.csv", sv_type=['TRS','NONTRS']), category="Step 2",
          subcategory="{model}",labels={"model": "{model}","figure": "some plot" })
    shell:
        """
        sudo mkdir -p {OUTPUT_DIR}/{sample_all}_ts/cyto_sv_ml &&
        {params.py39_dir} {main_dir}/software/CYTO-SV-ML/Pipeline_script/CYTO-SV-Auto-ML.py -s {sample_all} -o {OUTPUT_DIR}/{sample_all}_ts -k 5
        """             
