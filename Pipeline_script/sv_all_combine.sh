#!/bin/bash
main_dir=$1
sample=$2
size_k=$3

echo ${sample}
echo "# combine the SV annotation, sv vcf and svtyper info, sv complexity, and sv caller supp info" && date   
sed 's%_%:%g' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping | sed 's%sv:id%sv_id%g' > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp
awk 'FNR==NR{a[$1];b[$1]=$0;next} ($1 in a) {print $0"\t"b[$1]}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_anno
awk 'FNR==NR{a[$2];b[$2]=$0;next} ($2 in a) {print $0"\t"b[$2]}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_anno > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info 
awk 'FNR==NR{a[$6];b[$6]=$0;next} { if ($2 in a){print $0"\t"b[$2]} else {print $0"\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN\tNAN"}}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.kz.index_complex ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info >  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info.all_complex
awk 'FNR==NR{a[$3];b[$3]=$7;next} {if (FNR==1) {print $0"\tSUPP"} else if ($1 in a) {print $0"\t"b[$1]} else{print $0"\t1"}}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.sv_info.sim   ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info.all_complex > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.all_anno.all_info.all_complex.supp     
