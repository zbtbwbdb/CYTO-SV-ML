#!/bin/bash
main_dir=$1
cyto_dv_ml_dir=$2
sample=$3
size=$4
size_k=$((size/1000))

echo "# prepare the SV vcf files with size restriction and extract sv_vcf_info" && date 

# SV vcf transformation
for svtype in  breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta       
    do
        python ${cyto_dv_ml_dir}/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf f
        mv ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.re_id ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf    
        echo ${svtype} "SV size tf" && date     
        python ${cyto_dv_ml_dir}/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf $size down > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k
         echo ${svtype} "SV type tf" && date                
         python ${cyto_dv_ml_dir}/Pipeline_script/sv_vcf_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k
        echo ${svtype} "SV info tf" && date                 
        python ${cyto_dv_ml_dir}/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k
    done

cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf.${size_k}.trs_tf ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all
