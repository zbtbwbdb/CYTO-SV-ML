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
# all_callers_svtyper=['manta', 'delly', 'cnvnator', 'breakdancer']
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
        
# Run parliment2_sv
rule parliment2_sv:
#    singularity: 
#         "docker://docker.io/dongwonlee/parliament2-sing:v0.12"        
    input:         
        sample_bam=expand(INPUT_DIR+"/{sample}.bam",sample=SAMPLES)  
    output:  
        sample_vcf = expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf", sample=SAMPLES, sv_caller= parliment2_sv_callers)  
    threads: 8
    params:
        parliment_docker = parliment_docker,
        sm = SAMPLES
    shell:   
         """   
         bash {CYTO_SV_ML_DIR}/Pipeline_script/run_parliment2.sh {MAIN_DIR} {params.parliment_docker} {params.sm} 
         """    
        
# #  SV VCF preparation
rule sv_vcf_tf:
    input:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf", sample=SAMPLES, sv_caller=all_callers)
    output:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf", sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K, sv_type=['trs','nontrs']),        
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.sv_info.sim", sample=SAMPLES, size_k=SIZE_K, sv_caller=all_callers)
    params:
        sm = SAMPLES,   
    shell:
        """        
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_vcf_tf.sh {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {size}     
        """
        
# run sv merge
rule svmerge_qc:
    input:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf", sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K, sv_type=['trs','nontrs'])
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.{sv_type}_tf.all", sample=SAMPLES, size_k=SIZE_K, sv_type=['trs','nontrs']),           
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.{sv_type}", sample=SAMPLES, size_k=SIZE_K, sv_type=['trs','nontrs']),    
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all", size_k=SIZE_K, sample=SAMPLES)        
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {CYTO_SV_ML_DIR}/Pipeline_script/svmerge_qc.sh {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {SIZE_K}       
        """
        
# run svtyper qc
rule svtyper_qc:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.{sv_type}_tf.all", sample=SAMPLES, size_k=SIZE_K, sv_type=['trs','nontrs']),
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf", sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K, sv_type=['trs','nontrs'])        
    output:
        expand(OUTPUT_DIR+"/{sample}/${sample}.{size_k}k.{sv_type}_tf.all.svtyped.vcf.sv_info.sim", sample=SAMPLES, size_k=SIZE_K, sv_type=['trs','nontrs']),
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf.svtyped.vcf.sv_info.sim", sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K, sv_type=['trs','nontrs'])                       
    params:
        sm = SAMPLES,  
        sv_caller='@'.join(str(sc) for sc in all_callers)
    conda:
        "py27.yaml"          
    shell:
        """        
        bash {CYTO_SV_ML_DIR}/Pipeline_script/svtyper_qc.sh {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {params.sv_caller} {SIZE_K}     
        """
        
# run sv breakpoint sequence complexity       
rule sv_seq_complex:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all", sample=SAMPLES, size_k=SIZE_K)
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.bed.bpst_bpend.kz.index_complex", sample=SAMPLES, size_k=SIZE_K)
    params:
        sm = SAMPLES,
        py27_dir=config['py27_dir']      
    shell:
        """        
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_seq_complex.sh {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {params.py27_dir} {SIZE_K}     
        """

# run sv database annotation      
rule sv_database_ann:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.{sv_type}", sample=SAMPLES, size_k=SIZE_K, sv_type=['trs','nontrs'])   
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping.all_anno", sample=SAMPLES, size_k=SIZE_K)
    params:
        sm = SAMPLES,  
        py27_dir=config['py27_dir']
    shell:
        """        
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_database_ann.sh {MAIN_DIR} {CYTO_SV_ML_DIR} {params.sm} {params.py27_dir} {SIZE_K}      
        """

# run sv vcf info extraction          
rule sv_info_extract:
    input:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.sv_info.sim", sample=SAMPLES, size_k=SIZE_K, sv_caller=all_callers), 
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.{size_k}k.{sv_type}_tf.svtyped.vcf.sv_info.sim", sample=SAMPLES, sv_caller=all_callers, size_k=SIZE_K, sv_type=['trs','nontrs'])            
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping.all_info", sample=SAMPLES, size_k=SIZE_K)
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_info_extract.sh {MAIN_DIR} {params.sm} {SIZE_K}     
        """

# combine all sv features           
rule sv_all_combine:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.sv_id_mapping.{feature}", sample=SAMPLES, size_k=SIZE_K, feature=['all_anno','all_info']),
        expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.bed.bpst_bpend.kz.index_complex", sample=SAMPLES, size_k=SIZE_K)       
    output:
        report(expand(OUTPUT_DIR+"/{sample}/{sample}.{size_k}k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES, size_k=SIZE_K))  
    params:     
        sm = SAMPLES   
    shell:
        """     
        bash {CYTO_SV_ML_DIR}/Pipeline_script/sv_all_combine.sh {MAIN_DIR} {params.sm} {SIZE_K}  
        """
