#!/bin/bash
main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
size=$4
size_k=$((size/1000))

echo "# prepare the SV vcf files with size restriction and extract sv_vcf_info" && date 

# SV vcf transformation
for sv_caller in  breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta       
    do
        echo ${sv_caller}
        ls ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.list
        SURVIVOR merge ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.list 1000 0 1 0 0 10  ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s sc
        mv ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s.re_id ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s    
        echo ${sv_caller} "SV size tf" && date     
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s $size down > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k
        echo ${sv_caller} "SV type tf" && date                
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_vcf_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k
        echo ${sv_caller} "SV info tf" && date                 
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k
    done
