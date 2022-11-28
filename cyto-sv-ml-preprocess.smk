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
samples_information = pd.read_csv(config['sample_file'], sep='\t', header=None,index_col=False)
#print(samples_information)
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
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES)  
#        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf", sample=SAMPLES, sv_caller= all_callers)          
#         expand(OUTPUT_DIR+"/{sample}/{sample}.model_{analyses}.pdf", sample=SAMPLES, analyses=["onfusion_matrix.pdf", "aucroc_curve", "metrics"])
        
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
        sm = SAMPLES,  
        gd = GENDERS
    shell:
        """
         bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/run_chromoseq.sh {main_dir} {params.sm} {params.gd} 
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
        sm = SAMPLES
    shell:   
         """   
         bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/run_parliment2.sh {main_dir} {params.sm} 
         """    
        
# #  SV VCF preparation
rule sv_vcf_tf:
    input:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf", sample=SAMPLES, sv_caller=all_callers)
    output:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.10k.{sv_type}_tf", sample=SAMPLES, sv_caller=all_callers, sv_type=['trs','nontrs']),        
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.10k.sv_info.sim", sample=SAMPLES, sv_caller=all_callers)
    params:
        sm = SAMPLES,   
        size=size
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_vcf_tf.sh {main_dir} {params.sm} {params.size}     
        """
        
# run sv merge
rule svmerge_qc:
    input:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.10k.{sv_type}_tf", sample=SAMPLES, sv_caller=all_callers, sv_type=['trs','nontrs'])
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.{sv_type}", sample=SAMPLES, sv_type=['trs','nontrs']),
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all", sample=SAMPLES)        
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/svmerge_qc.sh {main_dir} {params.sm}     
        """
        
# run svtyper qc
rule svtyper_qc:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.{sv_type}_tf.all", sample=SAMPLES, sv_caller=all_callers_svtyper, sv_type=['trs','nontrs'])
    output:
        expand(OUTPUT_DIR+"/{sample}/${sample}.10k.{sv_type}_tf.all.svtyper.sv_info", sample=SAMPLES, sv_type=['trs','nontrs'])  
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/svtyper_qc.sh {main_dir} {params.sm}     
        """
        
# run sv breakpoint sequence complexity       
rule sv_seq_complex:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all", sample=SAMPLES)
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.bed.bpst_bpend.kz.index_complex", sample=SAMPLES)
    params:
        sm = SAMPLES   
    conda:
        "conda-py27.yaml"        
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_seq_complex.sh {main_dir} {params.sm}     
        """

# run sv database annotation      
rule sv_database_ann:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.{sv_type}", sample=SAMPLES, sv_type=['trs','nontrs'])   
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.sv_id_mapping.all_anno", sample=SAMPLES)
#     conda:
#         "conda-py27.yaml"
    params:
        sm = SAMPLES,  
        py27_dir=config['py27_dir']
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_database_ann.sh {main_dir} {params.sm} {params.py27_dir}      
        """

# run sv vcf info extraction          
rule sv_info_extract:
    input:
        expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf.10k.sv_info.sim", sample=SAMPLES, sv_caller=all_callers)  
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.sv_id_mapping.all_info", sample=SAMPLES)
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_info_extract.sh {main_dir} {params.sm}     
        """

# combine all sv features           
rule sv_all_combine:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.sv_id_mapping.{feature}", sample=SAMPLES,feature=['all_anno','all_info']),
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.bed.bpst_bpend.kz.index_complex", sample=SAMPLES)       
    output:
        report(expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_anno.all_info.all_complex.supp", sample=SAMPLES))  
    params:     
        sm = SAMPLES   
    shell:
        """     
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_all_combine.sh {main_dir} {params.sm}  
        """
            
    
