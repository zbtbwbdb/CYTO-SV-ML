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
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf c        
        ls ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.list        
        SURVIVOR merge ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.list 1000 0 1 0 0 10  ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s
        awk '($1!~"#"){split($3,b,":");print b[1]":"b[2]":"b[3]":"b[4]":"b[5]"\t"$3}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s.id
        awk 'FNR==NR{a[$1];c[$1]=$2;next}{split($3,b,":"); e=b[1]":"b[2]":"b[3]":"b[4]":"b[5]; if (($1!~"#")&&(e in a)) {$3=c[e]; print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.s.id ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id | sed 's% %\t%g' >   ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.s
        awk '($1~"#"){print $0}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id  > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.hd
        cat ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.hd ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id.s > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id
        echo ${sv_caller} "SV size tf" && date
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.re_id $size down > ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k
        echo ${sv_caller} "SV type tf" && date                
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_vcf_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k
        echo ${sv_caller} "SV info tf" && date                 
        python ${cyto_sv_ml_dir}/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k
    done
