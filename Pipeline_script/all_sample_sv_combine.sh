#!/bin/bash
main_dir=$1
sample_all=$2
SAMPLES_vector=$3

echo "# combine all sample sv " && date
for sample in $(echo ${SAMPLES_vector} | sed 's%@%%\t%g')
  do
    cat ${sample}.sv_all_combine} >> ${main_dir}/out/${sample_all}.sv.all.combine_all.tmp 
  done
awk '(FNR==1)||($1!~"sv_id"){print $0}' ${main_dir}/out/${sample_all}.sv.all.combine_all.tmp > ${main_dir}/out/${sample_all}.sv.all.combine_all
rm -rf  ${main_dir}/out/${sample_all}.sv.all.combine_all.tmp 
