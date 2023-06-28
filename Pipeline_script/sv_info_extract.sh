#!/bin/bash
main_dir=$1
sample=$2   
sv_caller_vector=$3
size_k=$4

echo "# starting SV info extraction" && date   
cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0
echo ${sv_caller_vector} | sed "s%@%\n%g" > ${main_dir}/out/${sample}/sv_caller_vector.tmp2
sc_ln=$(wc -l ${main_dir}/out/${sample}/sv_caller_vector.tmp2 | awk '{print $1}')

# # combine the sv_caller annotation vcf 
# n=0
# for i in $(seq 1 $sc_ln)
#    do
#        sv_caller=$(awk -v a="$i" '(FNR==a){print $1}' ${main_dir}/out/${sample}/sv_caller_vector.tmp2)  
#           echo $sv_caller "ok"      
#       if [ "${n}" == 0 ]; then
#           awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1
#           n=1
#       else
#           awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1    
#       fi
#       cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1 ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 
#    done
# awk '{print NF}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 | awk '!a[$1]++'   

# combine the svtyper annotation vcf 
n=0
for i in $(seq 1 $sc_ln)
   do
       sv_caller=$(awk -v a="$i" '(FNR==a){print $1}' sv_caller_vector.tmp)  
      if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim ]; then             
          echo ${sv_caller} "svtype ok" && date   
          if [ "${n}" == 0 ]; then                
              awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
              n=1                        
          else                    
              awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
          fi
      else   
          echo ${sv_caller} "svtype pass" && date     
          if [ -s ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.sv_info.sim ]; then                     
              cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.vcf.${size_k}k.sv_info.sim ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.all.svtyped.vcf.sv_info.sim.tmp
              if [ "${n}" == 0 ]; then
                  awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.svtyped.vcf.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
                  n=1
              else
                  awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.${sv_caller}.${size_k}k.svtyped.vcf.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t 
              fi   
          else
              echo ${sv_caller} "all pass" && date   
          fi
      fi
      cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp1t ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t          
done
    
echo "#check the complete status of sv_info extraction" && date  
awk '{print NF}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t | awk '!a[$1]++'   
info_check=$(awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t | awk '!a[$1]++' | wc -l )
if (( $info_check>1 )); then 
    echo "sv_info missing !!!"
    awk '{if (NF==2){print $0}}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info_missing       
    awk '{if (NF==2){print $0"\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t."} else{print $0}}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info   
    info_check2=$(awk '{print NF}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info | awk '!a[$1]++' | wc -l )
    if (( $info_check2>1 )); then   
        awk '{print NF}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info  | awk '!a[$1]++'       
        echo "sv_info2 missing !!!"    
        exit
    fi
else
    cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.tmp0t ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.sv_id_mapping.all_info             
fi
rm ${main_dir}/out/${sample}/${sample}.*tmp*      
rm ${main_dir}/out/${sample}/sv_caller_vector.tmp2
