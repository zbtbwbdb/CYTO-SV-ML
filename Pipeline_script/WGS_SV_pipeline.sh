#################################################################################################################
# pre-installation
# Parliment2: docker pull dongwonlee/parliament2-sing:v0.12
# ChromoSeq: docker pull zatawada/docker-basespace_chromoseq_v2:master
# SURVIVOR: https://github.com/fritzsedlazeck/SURVIVOR
# SVTyper: https://github.com/hall-lab/svtyper
# SeqComplex: https://github.com/caballero/SeqComplex
# Komplexity: https://github.com/eclarke/komplexity
#################################################################################################################

#!/bin/bash
main_dir=$1
sample_id_list=$2
sl=$(wc -l ${main_dir}/${sample_id_list} | awk '{print $1}')

# SV preprocess for each sample
for i in 1 #$(seq 1 $sl )
    do

        sample=$(awk -v a="$i" 'FNR==a{print $1}' ${main_dir}/${sample_id_list})
        gender=$(awk -v a="$i" 'FNR==a{print $2}' ${main_dir}/${sample_id_list})
        echo $i $sample $gender
# /usr/local/bin/aws s3 cp s3://gdr-prod-725905126021-app-mds-dr/reanalyzed_crams/${sample}.cram ${main_dir}/out/${sample}/

# /usr/local/bin/aws s3 cp s3://gdr-prod-725905126021-app-mds-dr/reanalyzed_crams/${sample}.cram.crai ${main_dir}/out/${sample}/

# # /usr/local/bin/aws s3 cp --recursive s3://gdr-prod-725905126021-app-mds-dr/results/reanalyzed/${sample}/svtyped_vcfs/ ${main_dir}/out/${sample}/

# sudo chmod 777 -R ${main_dir}/out/${sample}/*

# samtools view -b -@ 4 -m 2G -T ${main_dir}/reference/hg38/hs38.fasta ${main_dir}/out/${sample}/${sample}.cram -o ${main_dir}/out/${sample}/${sample}.bam

# rm -rf ${main_dir}/out/${sample}/${sample}.cram*

# samtools index -@ 4 -b ${main_dir}/out/${sample}/${sample}.bam

#         # SV Parliment2 run
#         sudo docker run --rm --privileged -v ${main_dir}/in/:/home/dnanexus/in -v ${main_dir}/out/${sample}/:/home/dnanexus/out docker.io/dongwonlee/parliament2-sing:v0.12 --bam ${sample}.bam  --bai ${sample}.bam.bai -r reference/hg38/hs38.fasta --fai reference/hg38/hs38.fasta.fai --prefix ${sample} --filter_short_contigs --lumpy
#         sudo docker run --rm --privileged -v ${main_dir}/in/:/home/dnanexus/in -v ${main_dir}/out/${sample}/:/home/dnanexus/out docker.io/dongwonlee/parliament2-sing:v0.12 --bam ${sample}.bam  --bai ${sample}.bam.bai -r reference/hg38/hs38.fasta --fai reference/hg38/hs38.fasta.fai --prefix ${sample} --filter_short_contigs --breakdancer --cnvnator --delly_deletion --delly_insertion --delly_inversion --delly_duplication

#        # SV ChromoSeq run
#         sed "s%XXXXXX%${sample}%g" ${main_dir}/software/docker-basespace_chromoseq/lsf/inputs.json | sed "s%Male%${gender}%g" > ${main_dir}/software/docker-basespace_chromoseq/lsf/inputs.json.tmp    
#         sudo docker run --rm --privileged  -v ${main_dir}/:/scratch  --entrypoint /bin/sh zatawada/docker-basespace_chromoseq_v2:master -c '/usr/bin/java -Dconfig.file=/scratch/software/docker-basespace_chromoseq/lsf/application.new.conf -jar /opt/cromwell-36.jar run -t wdl -i /scratch/software/docker-basespace_chromoseq/lsf/inputs.json.tmp /scratch/software/docker-basespace_chromoseq/workflow_files/Chromoseq.v17.wdl'

#       # SV vcf preparation
#         sudo mkdir ${main_dir}/out/${sample}/
#         sudo chmod -R 777 ${main_dir}/out/
#         cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.*.vcf  ${main_dir}/out/${sample}/
#         gunzip -f  ${main_dir}/out/${sample}/${sample}.svs_annotated.vcf.gz 
        ##awk '{if ($8~"SVTYPE=BND") {$9="GT:SU:PE:SR"; split($10,a,":"); if (length(a)==2) {a[3]="0,0"; split(a[2],c,",");split(a[3],d,","); b=c[2]+d[2];$10=a[1]":"b":"c[2]":"d[2];$8=$8";SU="b";PE="c[2]";SR="d[2]; print $0} else {split(a[2],c,",");split(a[3],d,","); b=c[2]+d[2];$10=a[1]":"b":"c[2]":"d[2];$8=$8";SU="b";PE="c[2]";SR="d[2]; print $0}}}' ${main_dir}/out/${sample}/${sample}.svs_annotated.vcf  | sed 's% %\t%g' > ${main_dir}/out/${sample}/${sample}.svs_annotated.vcf.cr
#         cp ${main_dir}/out/${sample}/${sample}.svs_annotated.vcf ${main_dir}/out/${sample}/${sample}.manta.vcf        
#        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/ichnorcnv_tf.py ${main_dir}/out/${sample}/testing.segs.txt ${main_dir}/out/${sample}/${sample}.ichnorcnv.vcf
#        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/${sample}.ichnorcnv.vcf f
#        mv ${main_dir}/out/${sample}/${sample}.ichnorcnv.vcf.re_id ${main_dir}/out/${sample}/${sample}.ichnorcnv.vcf 
#        python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_id_tf.py ${main_dir}/out/${sample}/${sample}.breakdancer.vcf f
#        mv ${main_dir}/out/${sample}/${sample}.breakdancer.vcf.re_id ${main_dir}/out/${sample}/${sample}.breakdancer.vcf 
#        cp ${main_dir}/out/${sample}/test123.sim.tf.rh.svtyper.vcf ${main_dir}/out/${sample}/${sample}.manta.svtyped.vcf 
         
#        # SV vcf consoldation
#           for svtype in manta ichnorcnv delly.deletion delly.deletion delly.duplication delly.inversion cnvnator breakdancer       
#               do
#                  echo ${svtype}
                 
#                  python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/${sample}.${svtype}.vcf 10000 down > ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k
#                  python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_vcf_tf.py ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k
#                  python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k  
#             done
#           for svtype in manta delly cnvnator breakdancer       
#               do     
#                  echo ${svtype}              
#                  python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_size.py ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf 10000 down > ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k
#                  python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k  
#               done
              
#         ls ${main_dir}/out/${sample}/${sample}.*.vcf.10k.nontrs_tf | grep -v svtyper > ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.list
#         SURVIVOR merge ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all 
#         ls ${main_dir}/out/${sample}/${sample}.*.vcf.10k.*trs_tf | grep -v svtyper > ${main_dir}/out/${sample}/${sample}.10k.all.list
#         SURVIVOR merge ${main_dir}/out/${sample}/${sample}.10k.all.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/${sample}.10k.sv.all 
#           python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py  ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all SUPP
#           python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_consolidate_id_mapping.py ${main_dir}/out/${sample}/${sample}.10k.sv.all          
          
# #         # SV svtyper run           
#         svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all -B ${main_dir}/out/${sample}/${sample}.bam > ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all.svtyper
#         python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all.svtyper
#          svtyper-sso --core 8 --max_reads 100000 -i ${main_dir}/out/${sample}/${sample}.manta.vcf.10k.trs_tf -B ${main_dir}/out/${sample}/${sample}.bam > ${main_dir}/out/${sample}/${sample}.10k.trs_tf.svtyper
#         python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_info_tf_sim.py ${main_dir}/out/${sample}/${sample}.10k.trs_tf.svtyper 
# rm -rf ${main_dir}/out/${sample}/${sample}.bam*

#         # SV sequence complexity run 
#         # make bed file for SV breakpoints
#          python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_vcf_bed_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all
#         awk 'FNR==NR{a[$1];b[$1]=$2;next}{c=b[$1]-150 ; if (($2>=150)&&($2<=c)) {$2=$2-150; $3=$2+150; print $1"\t"$2"\t"$3"\t"$5"\t"$6} else if ($2>c) {$2=c-150;$3=c+150; print $1"\t"$2"\t"$3"\t"$5"\t"$6} else if ($2<150){$2=1;$3=300; print $1"\t"$2"\t"$3"\t"$5"\t"$6}}' ${main_dir}/reference/hg38_chromosome_size.txt ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed | sed 's% %\t%g' > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst
#         awk 'FNR==NR{a[$1];b[$1]=$2;next}{$1=$4; $4=$1; c=b[$1]-150 ; if ($3>=c) {$3=c+150;$2=c-150; print $1"\t"$2"\t"$3"\t"$5"\t"$6} else {$2=$3-150; $3=$3+150;  print $1"\t"$2"\t"$3"\t"$5"\t"$6}}' ${main_dir}/reference/hg38_chromosome_size.txt ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed | sed 's% %\t%g' > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpend

# # #         export PATH=${main_dir}/software/SeqComplex:$PATH
# # #         cd ${main_dir}/software/SeqComplex
# # #          for bp in bpst bpend
# # #              do
# # # #                 echo ${bp}
# # # #                 # make bed file for SV breakpoints        
# # # #                 bedtools getfasta -fi ${main_dir}/reference/hg38/hs38.fasta -bed ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp} -fo ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out
                
# # # #                 # remove lowcomplex line with "NNNN"+ 1 "A/T/C/G"
# # # #                 grep NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out -n > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.cr
# # # #                 cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.tmp
# # # #                 for i in $(awk '{print $1}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.cr | cut -d ":" -f 1)
# # # #                     do
# # # #                         echo $i
# # # #                         i1=$((i-1))
# # # #                         awk -v a="$i1" -v b="$i" '{if ((FNR==a)||(FNR==b)) {print $0"\tlc"} else {print $0}}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.tmp > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.lc
# # # #                         cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.lc ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.tmp
# # # #                     done
# # # #                 awk '($NF!="lc"){print $0}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.tmp > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.lc
# # # #                # rm ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.tmp   
                
# # # #                 # SV breakpoint complexity
# # # #                 perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.lc
# # # #                 kz --fasta < ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.kz
# # #                 awk 'FNR==NR{a[$1]; b[$1]=$0; next}{c=$1":"$2"-"$3; if (c in a) {print $0"\t"b[c]}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.kz ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp} > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.kz.index
# # #                 awk 'FNR==NR{a[$1]; b[$1]=$0; next}($6 in a) {print $0"\t"b[$6]}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.complex ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.kz.index > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.${bp}.fa.out.kz.index_complex
# # #              done
#         awk 'FNR==NR{a[$5];b[$5]=$0;next} ($5 in a){print $0"\t"b[$5]}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst.fa.out.kz.index_complex ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpend.fa.out.kz.index_complex  > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst_bpend.fa.out.kz.index_complex
#          awk 'FNR==NR{a[$5];b[$5]=$0;next} ($6 in a){print $0"\t"b[$6]}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst_bpend.fa.out.kz.index_complex  ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed | awk '{if (FNR==1) {print "CHROM\tPOS\tEND\tsv_chr2\tsvtype\tID\tsv_bp_end_CHROM\tsv_bp_end_POS\tsv_bp_end_END\tsvtype\tID\tsv_bp_end_id\tsv_bp_end_length\tsv_bp_end_cc0\tsv_bp_end_cc1\tsv_bp_end_id\tsv_bp_end_cc_v1\tsv_bp_end_cc_v2\tsv_bp_end_cc_v3\tsv_bp_end_cc_v4\tsv_bp_end_cc_v5\tsv_bp_end_cc_v6\tsv_bp_end_cc_v7\tsv_bp_end_cc_v8\tsv_bp_end_cc_v9\tsv_bp_end_cc_v10\tsv_bp_end_cc_v11\tsv_bp_end_cc_v12\tsv_bp_end_cc_v13\tsv_bp_end_cc_v14\tsv_bp_end_cc_v15\tsv_bp_end_cc_v16\tsv_bp_end_cc_v17\tsv_bp_end_cc_v18\tsv_bp_end_cc_v19\tsv_bp_end_cc_v20\tsv_bp_end_cc_v21\tsv_bp_end_cc_v22\tsv_bp_end_cc_v23\tsv_bp_end_cc_v24\tsv_bp_st_CHROM\tsv_bp_st_POS\tsv_bp_st_END\tsvtype\tID\tsv_bp_st_id\tsv_bp_st_length\tsv_bp_st_cc0\tsv_bp_st_cc1\tsv_bp_st_id\tsv_bp_st_cc_v1\tsv_bp_st_cc_v2\tsv_bp_st_cc_v3\tsv_bp_st_cc_v4\tsv_bp_st_cc_v5\tsv_bp_st_cc_v6\tsv_bp_st_cc_v7\tsv_bp_st_cc_v8\tsv_bp_st_cc_v9\tsv_bp_st_cc_v10\tsv_bp_st_cc_v11\tsv_bp_st_cc_v12\tsv_bp_st_cc_v13\tsv_bp_st_cc_v14\tsv_bp_st_cc_v15\tsv_bp_st_cc_v16\tsv_bp_st_cc_v17\tsv_bp_st_cc_v18\tsv_bp_st_cc_v19\tsv_bp_st_cc_v20\tsv_bp_st_cc_v21\tsv_bp_st_cc_v22\tsv_bp_st_cc_v23\tsv_bp_st_cc_v24\n"$0} else {print $0}}' > ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst_bpend.kz.index_complex    
          
#         # SV vcf simplified transformation
         python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_vcf_sim.py ${main_dir}/out/${sample}/${sample}.10k.sv.all
         awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs_anno
         awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $5"\t"$1"\t"$2"\t"$3"\t"$1"\t"$4}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs_anno
        
#        # SV database annotation label
       for SV_database_name in 1000_gall # gnomad_gall control_gall cytoatlas_s cosmic_s control_g gnomad_g gnomad_qc 1000_g 
           do
               echo ${SV_database_name}
               if [ -s ${main_dir}/SV_database/${SV_database_name}.nontrs.gz ]; then
                   python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_database_mapping.py -i ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs -t ${main_dir}/SV_database/${SV_database_name}.nontrs.gz -d 1000 -p 0.7 -o ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}  
               else
                   awk 'FNR!=1{print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs > ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs.${SV_database_name}
               fi
               if [ -s ${main_dir}/SV_database/${SV_database_name}.trs ]; then
                   python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_bnd_database_mapping.py ${main_dir}/SV_database/${SV_database_name}.trs ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000 
               else
                   awk 'FNR!=1{print $0"\tNAN"}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs >  ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs.${SV_database_name}_1000 
               fi
                   # SV database annotation consolidation/transformation
#                    python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs ${SV_database_name}                    
#                    python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/sv_db_tf.py ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs ${SV_database_name}_1000  
           done

        # combine the sv annotation vcf 
#           cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0
#           n=0
#           for svtype in breakdancer cnvnator delly.deletion delly.duplication delly.inversion ichnorcnv manta 
#               do
#                     echo $svtype "ok"      
#                 if [ "${n}" == 0 ]; then
#                     awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1
#                     n=1
#                 else
#                     awk 'FNR==NR{a[$3];b[$3]=$0;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1    
#                 fi
#                 cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1 ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 
#              done
#           awk '{print NF}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 | awk '!a[$1]++'          
          
#           cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0 ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t
#           n=0
#           for svtype in breakdancer cnvnator delly ichnorcnv manta 
#                 do
#                 if [ -s ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim ]; then             
#                     echo $svtype "svtype ok"
#                     if [ "${n}" == 0 ]; then                
#                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
#                         n=1                        
#                     else                    
#                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
#                     fi
#                 else
#                     echo $svtype "svtype mock"                
#                     cp ${main_dir}/out/${sample}/${sample}.${svtype}.vcf.10k.sv_info.sim ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim.tmp
#                     if [ "${n}" == 0 ]; then
#                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ($2 in a) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
#                         n=1
#                     else
#                         awk 'FNR==NR{a[$3];b[$3]=$6;next} {if ((FNR!=1)&&($2 in a)) {print $0"\t"b[$2]} else {print $0}}' ${main_dir}/out/${sample}/${sample}.${svtype}.svtyped.vcf.10k.sv_info.sim.tmp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t 
#                     fi                     
#                 fi
#                 cp ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp1t ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t          
#               done
#                 awk '{if (NF==24){print $0"\t."} else{print $0}}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.tmp0t > ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_info
#                 awk '{print NF}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_info | awk '!a[$1]++'                  
                
#                #rm ${main_dir}/out/${sample}/${sample}.*tmp*
               
#          # combine the sv annotation and complexity and svtyper info               
#           cat ${main_dir}/out/${sample}/${sample}.10k.sv.all.trs_anno ${main_dir}/out/${sample}/${sample}.10k.sv.all.nontrs_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno  
#           awk 'FNR==NR{a[$1];b[$1]=$2;next} ($1 in a) {print b[$1]"\t"$0}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.sv_id_mapping 
#           awk 'FNR==NR{a[$2];b[$2]=$0;next} ($1 in a) {print $0"\t"b[$1]}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.sv_id_mapping.all_info ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.sv_id_mapping > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info 
#           awk 'FNR==NR{a[$6];b[$6]=$0;next} ($1 in a){print $0"\t"b[$1]}'  ${main_dir}/out/${sample}/${sample}.10k.sv.all.bed.bpst_bpend.kz.index_complex ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info >  ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex
#            awk 'FNR==NR{a[$3];b[$3]=$7;next} {if ($1 in a) {print $0"\t"b[$1]} else{print $0"\t1"}}'  ${main_dir}/out/${sample}/${sample}.10k.nontrs_tf.all.sv_info.sim   ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex > ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex.supp          
#            awk '{print NF}' ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_anno.all_info.all_complex.supp | awk '!a[$1]++'         
            done
#         # combine the sample SV into cohort dataset
#         sample_all="cohort_name" # please create your own cohort name here
#         cat ${main_dir}/out/${sample}/${sample}.10k.sv.all.all_combine >> ${main_dir}/out/${sample_all}.sv.all.all_combine

# # SV AutoML run
# python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/AutoML.py ${main_dir}/out/${sample_all}.sv.all.all_combine
# # Demo: python AutoML.py example/input.csv