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
print(SAMPLES_vector)

MAIN_DIR = config['main_dir']
INPUT_DIR = config['main_dir']+'/in'
OUTPUT_DIR = config['main_dir']+'/out'
LOG_DIR = config['main_dir']+'/out/log'
CYTO_SV_ML_DIR = config['cyto_sv_ml_dir']
SOFTWARE_DIR = config['cyto_sv_ml_dir']+'/software'
DATABASE_DIR = config['cyto_sv_ml_dir']+'/SV_database'
CYTO_BAND_FILE = config['cyto_band_file']

# build R-shiny based user interface in docker container
# 2 GB space is required to build and run this docker image 
rule interface_docker:
    input:
        expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}_{sv_type}_{k}_ts_EXP", cohort_name=cohort_name, k=[0,1,2,3,4], sv_type=['trs','nontrs'])   
    output:
        expand(OUTPUT_DIR+"/{cohort_name}/cyto_sv_ml/{cohort_name}.all_pred", cohort_name=cohort_name),
        CYTO_SV_ML_DIR+"/docker/cyto-sv-ml/data/cyto_sv_ml.RData"
    params:
        kfolds=config['kfolds']         
    shell:
        """  
        echo {input} && bash {CYTO_SV_ML_DIR}/Pipeline_script/interface_docker.sh {MAIN_DIR} {CYTO_SV_ML_DIR} {cohort_name} {CYTO_BAND_FILE} {params.kfolds} 
        """      
