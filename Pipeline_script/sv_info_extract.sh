#!/bin/bash
main_dir=$1
sample=$2   

#        extract the sv vcf info columns 
           cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0
           n=0
           for svtype in breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta 
               do
                     echo $svtype "ok"      
                 if [ "${n}" == 0 ]; then
                     awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1
                     n=1
                 else
                     awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1    
                 fi
                 cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1 ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 
              done
           awk '{print NF}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 | awk '!a[$1]++'          
          
           cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t
           n=0
           for svtype in breakdancer cnvnator delly ichnorcnv manta 
                 do
                 if [ -s ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim ]; then             
                     echo $svtype "svtype ok"
                     if [ "${n}" == 0 ]; then                
                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
                         n=1                        
                     else                    
                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
                     fi
                 else
                     echo $svtype "svtype mock"                
                     cp ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim.tmp
                     if [ "${n}" == 0 ]; then
                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
                         n=1
                     else
                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
                     fi                     
                 fi
                 cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t          
               done
                 awk '{if (NF==24){print $0"\t."} else{print $0}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_info
                rm ${main_dir}/out/${sample}/${sample}.*tmp*                 
                 awk '{print NF}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_info | awk '!a[$1]++'                  
               