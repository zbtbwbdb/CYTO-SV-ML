#!/bin/bash
main_dir=$1
sample=$2   
py27_dir=$3

# intiate the touch of anno files
         awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs_anno
         awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $5"\t"$1"\t"$2"\t"$3"\t"$1"\t"$4}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs_anno
# SV database annotation label
       for SV_database_name in 1000_gall gnomad_gall control_gall cytoatlas_s cosmic_s control_g gnomad_g gnomad_qc 1000_g 
           do
               echo ${SV_database_name}
               if [ -s ${main_dir}/SV_database/${SV_database_name}.nontrs.gz ]; then
                   ~/mambaforge/envs/py27/bin/python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_database_mapping.py -i ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs -t ${main_dir}/SV_database/${SV_database_name}.nontrs.gz -d 1000 -p 0.7 -o ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}  
               else
                   awk 'FNR!=1{print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}
               fi
               if [ -s ${main_dir}/SV_database/${SV_database_name}.trs ]; then
                   python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_bnd_database_mapping.py ${main_dir}/SV_database/${SV_database_name}.trs ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000 
               else
                   awk 'FNR!=1{print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs >  ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs.${SV_database_name}_1000 
               fi
#                    SV database annotation consolidation/transformation
                    python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs ${SV_database_name}                    
                    python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000  
           done
           cat ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs_anno ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno  
           awk 'FNR==NR{a[$1];b[$1]=$2;next} ($1 in a) {print b[$1]"\t"$0}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_anno            