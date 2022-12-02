#!/bin/bash
main_dir=$1
sample=$2   
size_k=$3

echo "# starting SV info extraction" && date   
cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t
n=0
for sv_callcer in breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta 
  do
      if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.${size_k}k.all.svtyped.vcf.sv_info.sim ]; then             
          echo ${sv_callcer} "svtype ok" && date   
          if [ "${n}" == 0 ]; then                
              awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.${size_k}k.all.svtyped.vcf.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
              n=1                        
          else                    
              awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.${size_k}k.all.svtyped.vcf.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
          fi
      else   
          echo ${sv_callcer} "svtype pass" && date     
          if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.vcf.${size_k}k.sv_info.sim ]; then                     
              cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.vcf.${size_k}k.sv_info.sim ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.${size_k}k.all.svtyped.vcf.sv_info.sim.tmp
              if [ "${n}" == 0 ]; then
                  awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.${size_k}k.svtyped.vcf.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
                  n=1
              else
                  awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_callcer}.${size_k}k.svtyped.vcf.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
              fi   
          else
              echo ${sv_callcer} "all pass" && date   
          fi
      fi
      cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t          
done
    
echo "# consolidate SV info" && date    
awk '{if (NF==24){print $0"\t."} else{print $0}}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info
rm ${main_dir}/out/${sample}/${sample}.*tmp*                 
 
# check the complete status of sv_info extraction
awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info | awk '!a[$1]++'       
info_check=$(awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info | awk '!a[$1]++'  | wc -l | awk '{print $1}')
if (($info_check>1)); then echo "sv_info missing !!!"; fi    
