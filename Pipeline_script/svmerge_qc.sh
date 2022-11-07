#!/bin/bash
main_dir=$1
sample=$2

# SV vcf consoldation
ls ${main_dir}/out/${sample}/${sample}.*.vcf.10k.nontrs_tf | grep -v svtyper > ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.list
SURVIVOR merge ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all 
ls ${main_dir}/out/${sample}/${sample}.*.vcf.10k.*trs_tf | grep -v svtyper > ${main_dir}/out/${sample}/${sample}.10k.all.list
SURVIVOR merge ${main_dir}/out/${sample}/${sample}.10k.all.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/${sample}.10k.sv.all 
python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py  ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all SUPP
python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_consolidate_id_mapping.py ${main_dir}/out/${sample}/${sample}.10k.sv.all 

         # SV vcf simplified transformation
         python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_vcf_sim.py ${main_dir}/out/${sample}/${sample}.10k.sv.all
