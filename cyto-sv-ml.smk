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
print(samples_information)
samples_information.columns=['id','sex']
SAMPLES = list(samples_information['id'])
GENDERS = list(samples_information['sex'])

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
all_callers=parliment2_sv_callers+chromoseq_sv_callers
all_callers_svtyper=['manta', 'delly', 'cnvnator', 'breakdancer']

report: OUTPUT_DIR+"report/workflow.rst"
    
rule all:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.model_{analyses}.pdf", sample=SAMPLES, analyses=["onfusion_matrix.pdf", "aucroc_curve", "metrics"])

# Run parliment2_sv
rule parliment2_sv:
#    singularity: 
#         "docker://docker.io/dongwonlee/parliament2-sing:v0.12"        
    input:         
        sample_bam=expand(INPUT_DIR+"/{sample}.bam",sample=SAMPLES)  
    output:  
        sample_vcf = expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf",sample=SAMPLES , sv_caller= 'delly.duplication')  
    params:
        main_dir = config['main_dir'],
        sm = SAMPLES
    shell:   
         """   
         sudo singularity exec --bind {params.main_dir}/in:/home/dnanexus/in  --bind {params.main_dir}/out/{params.sm}:/home/dnanexus/out  --writable-tmpfs parliament2-sing_v0.12.sif /home/dnanexus/parliament2.py --bam {params.sm}.bam  --bai {params.sm}.bam.bai -r reference/hg38/hs38.fasta --fai reference/hg38/hs38.fasta.fai --prefix {params.sm} --filter_short_contigs --delly_duplication
         """    
        
# Run chromoseq_sv
rule chromoseq_sv:
#    singularity: 
#        "docker://docker.io/zatawada/docker-basespace_chromoseq_v2:master"
    input:
        sample_cram=expand(INPUT_DIR+"/{sample}.cram",sample=SAMPLES) 
    output:
        sample_vcf = expand(OUTPUT_DIR+"/{sample}/sv_caller_results/{sample}.{sv_caller}.vcf"), sv_caller=chromoseq_sv_callers)
    params:
        main_dir = config['main_dir'],
        sample = SAMPLES  
        gender = GENDERS
    shell:
        """
        sed "s%XXXXXX%{sample}%g" SOFTWARE_DIR/docker-basespace_chromoseq/lsf/inputs.json | sed "s%Male%${gender}%g" > SOFTWARE_DIR/docker-basespace_chromoseq/lsf/inputs.json.tmp,
        /usr/bin/java -Dconfig.file={SOFTWARE_DIR}/docker-basespace_chromoseq/lsf/application.new.conf -jar /opt/cromwell-36.jar run -t wdl -i SOFTWARE_DIR/docker-basespace_chromoseq/lsf/inputs.json.tmp SOFTWARE_DIR/docker-basespace_chromoseq/workflow_files/Chromoseq.v17.wdl     
        gunzip -f OUTPUT_DIR/{sample}/{sample}.svs_annotated.vcf.gz 
        cp OUTPUT_DIR/{sample}/{sample}.svs_annotated.vcf OUTPUT_DIR/{sample}/sv_caller_results/{sample}.manta.vcf    
        python SOFTWARE_DIR/CYTO-SV-ML/Pipeline_script/ichnorcnv_tf.py OUTPUT_DIR/{sample}/{sample}.segs.txt OUTPUT_DIR/{sample}/sv_caller_results/{sample}.ichnorcnv.vcf 
        """
        
# #  SV VCF preparation
rule sv_vcf_tf:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{sv_caller}.vcf", sample=SAMPLES, sv_caller=all_callers)
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{sv_caller}.vcf.10k.sv_info.sim", sample=SAMPLES, sv_caller=all_callers)
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_vcf_tf.sh {main_dir} {params.sm}     
        """
        
# run sv merge
rule svmerge_qc:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{sv_caller}.vcf.10k.{sv_type}_tf", sample=SAMPLES, sv_caller=all_callers, sv_type=['trs','nontrs'])
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.{sv_type}_anno", sample=SAMPLES, sv_type=['trs','nontrs'])
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/svmerge_qc.sh {main_dir} {params.sm}     
        """
        
# run svtyper qc
rule svtyper_qc:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{sv_caller}.svtyped.vcf", sample=SAMPLES, sv_caller=all_callers_svtyper)  
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{sv_caller}.svtyped.vcf.10k.sv_info", sample=SAMPLES, sv_caller=all_callers_svtyper)  
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
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_seq_complex.sh {main_dir} {params.sm}     
        """

# run sv database annotation      
rule sv_database_ann:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.{sv_type}_anno", sample=SAMPLES, sv_type=['trs','nontrs'])   
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.sv_id_mapping.all_anno", sample=SAMPLES)
#     conda:
#         "conda-py27.yaml"
    params:
        sm = SAMPLES,  
        py27_dir=config['py27_dir']
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_database_ann.sh {main_dir} {params.sm}     
        """

# run sv vcf info extraction          
rule sv_info_extract:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.{sv_caller}.vcf.10k.sv_info.sim", sample=SAMPLES, sv_caller=all_callers)  
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
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.bed.bpst_bpend.kz.index_complex", sample=SAMPLES),
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.nontrs_tf.all.sv_info.sim", sample=SAMPLES)        
    output:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_combine.supp", sample=SAMPLES)  
    params:
        sm = SAMPLES         
    shell:
        """        
        bash {SOFTWARE_DIR}/CYTO-SV-ML/Pipeline_script/sv_all_combine.sh {main_dir} {params.sm}     
        """
        
#  checkpoint for all sample sv data   
checkpoint all_sample_sv_ready:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_combine.supp", sample=SAMPLES)  
    output: "sample_sv_ready.out"
    run: 
        shell('echo {input} >> {OUTPUT_DIR}/log/sample_sv_ready.out')

def check_sample_file(*wildcards):
     return checkpoints.all_sample_sv_ready.get().output        

# combine all sample sv           
rule sv_info_combine:
    input:
        expand(OUTPUT_DIR+"/{sample}/{sample}.10k.sv.all.all_combine.supp", sample=SAMPLES)       
    output:
        expand(OUTPUT_DIR+"/{sample_all}.sv.all.combine_all", sample_all=sample_all)   
    conda:
        "conda_py27.yaml"
    params:
        sample_all = sample_all         
    shell:
        """        
         cat {main_dir}/out/*/*.10k.sv.all.all_combine >> {main_dir}/out/{sample_all}.sv.all.all_combine  
        """               
               
# # # run cyto-sv-ml model           
# # rule cyto_sv_ml:
#     output:
#         report(
#             expand(OUTPUT_DIR+"/{sample}/{sample}.model_{analyses}.pdf", sample=SAMPLES, analyses=["type_class_summary",score_metrics", "confusion_matrix", "aucroc_curve"])
#                   )
#     run:
#         python {main_dir}/software/CYTO-SV-ML/Pipeline_script/AutoML.py {main_dir}/out/{sample_all}.sv.all.all_combine      
    