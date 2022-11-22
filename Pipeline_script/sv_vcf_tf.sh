#!/bin/bash
main_dir=$1
sample=$2
size=$3

echo "# prepare the SV vcf files with size restriction and extract sv_vcf_info" && date 

# SV vcf transformation
for svtype in  breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta       
    do
        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf f
        mv ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.re_id ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf    
        echo ${svtype} "SV size tf" && date     
        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf $size down > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.10k
         echo ${svtype} "SV type tf" && date                
         python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_vcf_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.10k
        echo ${svtype} "SV info tf" && date                 
        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.10k
    done

cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf.10k.trs_tf ${main_dir}/out/${sample}/${sample}.10k.trs_tf.all
