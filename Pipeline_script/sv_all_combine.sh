#!/bin/bash
main_dir=$1
sample=$2

echo ${sample}
echo "# combine the SV annotation, complexity, vcf and svtyper info" && date               
awk 'FNR==NR{a[$2];b[$2]=$0;next} ($1 in a) {print $0"\t"b[$1]}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_info ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info 
awk 'FNR==NR{a[$6];b[$6]=$0;next} ($1 in a){print $0"\t"b[$1]}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst_bpend.kz.index_complex ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info >  ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex
awk 'FNR==NR{a[$3];b[$3]=$7;next} {if ($1 in a) {print $0"\t"b[$1]} else{print $0"\t1"}}'  ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all.sv_info.sim   ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex.supp          
