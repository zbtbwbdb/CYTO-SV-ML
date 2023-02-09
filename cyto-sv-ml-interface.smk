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

# build R-shiny based user interface in docker container
# 2 GB space is required to build and run this docker image 
rule interface_docker:
    input:
        OUTPUT_DIR+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_trs_EXP,
        OUTPUT_DIR+'/cyto_sv_ml/'+str(cohort_name)+'_'+sv_type+'_'+str(k)+'_nontrs_EXP
    output:
        expand(OUTPUT_DIR+"/{cohort_name}/{cohort_name}.Dockfile", cohort_name=cohort_name)        
    shell:
        """  
        cat {input} && bash {CYTO_SV_ML_DIR}/Pipeline_script/interface_docker.sh {MAIN_DIR} {cohort_name} {input} 
        """      
