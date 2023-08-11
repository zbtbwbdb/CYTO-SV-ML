#!/bin/bash
main_dir=$1
cyto_sv_ml_dir=$2
sample=$3
sv_caller_vector=$4
size_k=$5
svtyper="svtyper-sso --core 8"
#svtyper=${main_dir}"/mambaforge/envs/py27/bin/svtyper-sso --core 8"

echo ${sv_caller_vector} | sed "s%@%\n%g" > sv_caller_vector.tmp
sc_ln=$(wc -l sv_caller_vector.tmp | awk '{print $1}')
echo "# run sytyper for all sv callers" && date
for i in $(seq 1 $sc_ln)
   do
       sv_caller=$(awk -v a="$i" '(FNR==a){print $1}' sv_caller_vector.tmp)  
       echo ${sv_caller}   
       python ${cyto_sv_ml_dir}/Pipeline_script/trs_svtyper_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.trs_tf 
       ${svtyper} --max_reads 100000 -i ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.trs_tf.tmp -B ${main_dir}/in/${sample}.bam > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf      
       ${svtyper} --max_reads 100000 -i ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.nontrs_tf -B ${main_dir}/in/${sample}.bam > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.nontrs_tf.svtyped.vcf     
       if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf ]; then    
            awk '{if ($1~"#") {print $0} else if ($5~"\\["){split($5,a,"["); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0} else if ($5~"\\]"){split($5,a,"]"); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.re_id
            awk '($1!~"#"){if ($5~"\\["){split($5,a,"["); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0} else if ($5~"\\]"){split($5,a,"]"); split(a[2],b,":"); $3=$1":"$2":"b[2]":"b[1]":BND:"$3; print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf | sed 's% %\t%g' > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp                
       fi
       if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp ]; then      
            cat ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.nontrs_tf.svtyped.vcf ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf
            sudo rm -rf ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.trs_tf.svtyped.vcf.tmp
       else
            cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.nontrs_tf.svtyped.vcf ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf 
       fi
       python ${cyto_sv_ml_dir}/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf c
       sudo mv ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.re_id ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf
       python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf   
    done
sudo rm -rf sv_caller_vector.tmp 

# # SV svtyper run    
# # echo "# run sytyper for all nontrs SV" && date
# ${svtyper} --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all -B ${main_dir}/in/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.svtyped.vcf

# # echo "# run sytyper for all trs SV" && date
# python ${cyto_sv_ml_dir}/Pipeline_script/trs_svtyper_tf.py ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all
# ${svtyper} --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.tmp -B ${main_dir}/in/${sample}.bam > ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf

# awk '($1!~"#"){print $0}' ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf > ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf.tmp 
# cat ${main_dir}/out/${sample}/${sample}.${size_k}k.nontrs_tf.all.svtyped.vcf ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf.tmp > ${main_dir}/out/${sample}/${sample}.${size_k}k.all.svtyped.vcf
# sudo rm -rf ${main_dir}/out/${sample}/${sample}.${size_k}k.trs_tf.all.svtyped.vcf.tmp
# python ${cyto_sv_ml_dir}/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/${sample}.${size_k}k.all.svtyped.vcf c
# sudo mv ${main_dir}/out/${sample}/${sample}.${size_k}k.all.svtyped.vcf.re_id ${main_dir}/out/${sample}/${sample}.${size_k}k.all.svtyped.vcf
# python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${size_k}k.all.svtyped.vcf

sudo rm -rf ${main_dir}/in/${sample}.bam*
