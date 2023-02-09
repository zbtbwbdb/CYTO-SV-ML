import os
import sys
import pandas as pd
import numpy as np
import pathlib
import snakemake.io
from snakemake.utils import validate
from typing import Dict, Union, List

configfile: "config.yaml"

SAMPLES = config['sample']
GENDERS = config['gender']

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
SV_DB=config['sv_db']
size=int(config['size'])
SIZE_K=round(size/1000)
#report: OUTPUT_DIR+"/report/workflow.rst"
    
rule all:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES, size_k=SIZE_K)  
        
# Run chromoseq_sv
rule chromoseq_sv:
#    singularity: 
#        "docker://docker.io/zatawada/docker-basespace_chromoseq_v2:master"
    input:
        sample_cram=expand(INPUT_DIR+"/{sample}.cram",sample=SAMPLES) 
    output:
        sample_vcf = expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf", sample=SAMPLES, sv_caller=chromoseq_sv_callers)
    threads: 8        
    params:
        chromoseq_docker = chromoseq_docker,        
        sm = SAMPLES,  
        gd = GENDERS
    shell:
        """
         bash {CYTO_SV_ML_DIR}/Pipeline_script/run_chromoseq.sh {MAIN_DIR} {CYTO_SV_ML_DIR} {params.chromoseq_docker} {params.sm} {params.gd} 
        """
        
