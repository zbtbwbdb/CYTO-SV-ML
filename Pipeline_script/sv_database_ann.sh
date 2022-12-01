#!/bin/bash
main_dir=$1
cyto_dv_ml_dir=$2
sample=$3   
py27_dir=$4

echo "# intiate the touch of SV anno files" && date
awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs_anno
awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $5"\t"$1"\t"$2"\t"$3"\t"$1"\t"$4}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs_anno

echo "# SV database annotation label" && date
# SV uwstl sv label
for SV_database_name in uwstl_s 
    do
        echo ${SV_database_name} "ok" && date            
        if [ -s ${main_dir}/out/${sample}/${sample}.${SV_database_name}.nontrs.gz ]; then                  
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_database_mapping.py -i ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs -t ${main_dir}/out/${sample}/${sample}.${SV_database_name}.nontrs.gz -d 1000 -p 0.7 -o ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}  
        else
            awk 'FNR!=1{$4=$1"\t"$4; print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs | sed 's% %\t%g' > ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}
        fi
        if [ -s ${main_dir}/out/${sample}/${sample}/${SV_database_name}.trs ]; then
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_bnd_database_mapping.py ${main_dir}/out/${sample}./${sample}.${SV_database_name}.trs ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000 
        else
            awk 'FNR!=1{print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs | sed 's% %\t%g' >  ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs.${SV_database_name}_1000 
        fi
#           SV database annotation consolidation/transformation
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs ${SV_database_name} f                   
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000 f 
   done
   
# SV database annotation label
for SV_database_name in  1000_g 1000_gall control_g control_gall cosmic_s cytoatlas_s gnomad_g2ab gnomad_gall gnomad_g gnomad_qc centromere_qc dgv_g
    do
        echo ${SV_database_name} "ok" && date  
        if [ -s ${cyto_dv_ml_dir}/SV_database/${SV_database_name}.nontrs.gz ]; then
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_database_mapping.py -i ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs -t ${cyto_dv_ml_dir}/SV_database/${SV_database_name}.nontrs.gz -d 1000 -p 0.7 -o ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}  
        else
            awk 'FNR!=1{$4=$1"\t"$4; print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs | sed 's% %\t%g' > ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}
        fi
        if [ -s ${cyto_dv_ml_dir}/SV_database/${SV_database_name}.bp.trs ]; then
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_bnd_database_mapping.bp.py ${cyto_dv_ml_dir}/SV_database/${SV_database_name}.bp.trs ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000 
        else
            awk 'FNR!=1{print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs | sed 's% %\t%g' >  ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs.${SV_database_name}_1000 
        fi
           # SV database annotation consolidation/transformation
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs ${SV_database_name} f                   
            ${py27_dir} ${cyto_dv_ml_dir}/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000 f 
    done           

echo "# prepare SV DB annotation file" && date
cat ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs_anno ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno 
awk 'FNR==NR{a[$1];b[$1]=$2;next} ($1 in a) {print b[$1]"\t"$0}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_anno            
