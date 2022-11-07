#!/bin/bash
main_dir=$1
sample=$2
           for svtype in manta delly cnvnator breakdancer       
               do     
                  echo ${svtype}              
                  python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf 10000 down > ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k
                  python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k  
               done
               
 #         # SV svtyper run           
         # svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all -B ${main_dir}/out/${sample}/${sample}.bam > ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all.svtyper
         #python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all.svtyper
         # svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.manta.vcf.10k.trs_tf -B ${main_dir}/out/${sample}/${sample}.bam > ${main_dir}/out/${sample}/${sample}.10k.trs_tf.svtyper
         #python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.10k.trs_tf.svtyper 
 rm -rf ${main_dir}/out/${sample}/${sample}.bam*
