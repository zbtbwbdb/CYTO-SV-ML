#!/bin/bash
main_dir=$1
sample_all=$2
SAMPLES_vector=$3

echo "# combine all sample sv " && date
sudo rm -rf ${main_dir}/out/${sample_all}/
sudo mkdir ${main_dir}/out/${sample_all}/
sudo chmod 777 -R ${main_dir}/out/${sample_all}/
sed "s%@%\n%g" ${SAMPLES_vector} > ${SAMPLES_vector}.tmp
sm_ln=$(wc -l ${SAMPLES_vector}.tmp | awk '{print $1}')

for i in $(seq 1 $sm_ln)
  do
    sample=$(awk -v a="$i" '(FNR==a){print $1}' ${SAMPLES_vector}.tmp)
    sample_file="${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex.supp"
    echo ${sample} ${sample_file}
    if (( i==1 )); then
        awk -v a="$sample" '{if (FNR==1) {print $0"\tsample_id"} else {print $0"\t"a}}' ${sample_file} > ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all
    else
        awk -v a="$sample" '(FNR>1){print $0"\t"a}' ${sample_file} >> ${main_dir}/out/${sample_all}/${sample_all}.sv.all.combine_all      
        fi
  done
sudo rm -rf ${SAMPLES_vector}.tmp

