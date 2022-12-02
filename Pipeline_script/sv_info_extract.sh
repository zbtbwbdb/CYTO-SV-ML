#!/bin/bash
main_dir=$1
sample=$2   
size_k=$3

cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0
n=0
for svtype in breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta 
    do
      if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k.sv_info.sim ]; then 
          echo $svtype "ok" && date       
          if [ "${n}" == 0 ]; then
              awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1
              n=1
          else
              awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1    
          fi
        else
            echo $svtype "pass" && date       
        fi
        cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1 ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0 
       done

# check the complete status of sv_info extraction
awk '{print NF}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0 | awk '!a[$1]++'          
info_check=$(awk '{print NF}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0 | awk '!a[$1]++' | wc -l | awk '{print $1}')
if (($info_check>1)); then echo "sv_info missing !!!"; fi

cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0 ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t
n=0
for svtype in breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta 
  do
      if [ -s ${main_dir}/out/${sample}/${sample}.${svtype}.${size_k}k.svtyped.vcf.sv_info.sim ]; then             
          echo $svtype "svtype ok" && date   
          if [ "${n}" == 0 ]; then                
              awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.${size_k}k.svtyped.vcf.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
              n=1                        
          else                    
              awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.${size_k}k.svtyped.vcf.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
          fi
      else   
          echo $svtype "svtype pass" && date     
          if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k.sv_info.sim ]; then                     
              cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.${svtype}.vcf.${size_k}k.sv_info.sim ${main_dir}/out/${sample}/${sample}.${svtype}.${size_k}k.svtyped.vcf.sv_info.sim.tmp
              if [ "${n}" == 0 ]; then
                  awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.${size_k}k.svtyped.vcf.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
                  n=1
              else
                  awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.${size_k}k.svtyped.vcf.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
              fi   
          else
              echo $svtype "all pass" && date   
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
