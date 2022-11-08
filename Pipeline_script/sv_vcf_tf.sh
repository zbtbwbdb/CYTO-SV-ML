#!/bin/bash
main_dir=$1
sample=$2

echo "# prepare the SV vcf files (manta + ichnorcnv)" && date 
gunzip -f ${main_dir}/out/${sample}/sv_caller_results/${sample}.svs_annotated.vcf.gz 
cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.svs_annotated.vcf ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf    
python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/ichnorcnv_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.segs.txt ${main_dir}/out/${sample}/sv_caller_results/${sample}.ichnorcnv.vcf 

# SV vcf transformation
for svtype in  breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta       
    do
        echo ${svtype} "SV size tf" && date
        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf 10000 down > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.10k
        echo ${svtype} "SV type tf" && date                
        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_vcf_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.10k
        echo ${svtype} "SV info tf" && date                 
        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.10k  
    done
cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf.10k.trs_tf ${main_dir}/out/${sample}/${sample}.10k.trs_tf.all
