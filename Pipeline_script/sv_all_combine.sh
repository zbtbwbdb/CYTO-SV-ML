#!/bin/bash
main_dir=$1
sample=$2

echo ${sample}
echo "# combine the SV annotation, complexity, vcf and svtyper info" && date               
awk 'FNR==NR{a[$1];b[$1]=$0;next} ($1 in a) {print $0"\t"b[$1]}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_anno
awk 'FNR==NR{a[$1];b[$1]=$0;next} ($1 in a) {print $0"\t"b[$1]}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_info ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info 
awk 'FNR==NR{a[$6];b[$6]=$0;next} { if ($2 in a){print $0"\t"b[$2]} else {print $0"\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN"}}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst_bpend.kz.index_complex ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info >  ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex
awk 'FNR==NR{a[$3];b[$3]=$7;next} {if (FNR==1) {print $0"\tSUPP"} else if ($1 in a) {print $0"\t"b[$1]} else{print $0"\t1"}}'  ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all.sv_info.sim   ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex.supp     
