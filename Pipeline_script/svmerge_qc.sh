#!/bin/bash
main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
size_k=$4

echo "# merge all non/trs SV" && date
ls ${main_dir}/out/${sample}/sv_caller_results/${sample}.*.vcf.${size_k}k.nontrs_tf | grep -v svtyper > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.nontrs_tf.list
SURVIVOR merge ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.nontrs_tf.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all
cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf.${size_k}k.trs_tf ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all

echo "# merge all SV" && date
ls ${main_dir}/out/${sample}/sv_caller_results/${sample}.*.vcf.${size_k}k.*trs_tf | grep -v svtyper > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.all.list
SURVIVOR merge ${main_dir}/out/${sample}/sv_caller_results/${sample}.${size_k}k.all.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all 

echo "# SV caller SUPP info extraction" && date
python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py  ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all SUPP
python ${cyto_sv_ml_dir}/Pipeline_script/sv_consolidate_id_mapping.py ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all 

echo "# SV vcf simplified transformation" && date
python ${cyto_sv_ml_dir}/Pipeline_script/sv_vcf_sim.py ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all
