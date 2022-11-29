#!/bin/bash
main_dir=$1
sample=$2   
gender=$3
echo $sample $gender
echo "# set up the configure file for chromoseq" && date
sed "s%/scratch/out/XXXXXX%/scratch/out/${sample}/sv_caller_results/%g" ${main_dir}/software/docker-basespace_chromoseq/lsf/inputs.json | sed "s%XXXXXX%${sample}%g"  | sed "s%Male%${gender}%g" > ${main_dir}/software/docker-basespace_chromoseq/lsf/inputs.json.tmp

echo "# run chromoseq docker" && date
 
sudo docker run --rm --privileged  -v ${main_dir}/:/scratch  --entrypoint /bin/sh  docker.io/zatawada/docker-basespace_chromoseq_v2:master -c '/usr/bin/java -Dconfig.file=/scratch/software/docker-basespace_chromoseq/lsf/application.new.conf -jar /opt/cromwell-36.jar run -t wdl -i /scratch/software/docker-basespace_chromoseq/lsf/inputs.json.tmp /scratch/software/docker-basespace_chromoseq/workflow_files/Chromoseq.v17.wdl'     

echo "# prepare the SV db uwstl_s file for chromoseq pipeline" && date 
awk '($1=="DUP")||($1=="DEL")||($1=="BND")||($1=="INV"){print $2"\t"$3"\t"$5"\t"$1"\t"$2":"$3":"$5":"$1}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.chromoseq.txt > ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.nontrs
sort -k 1,1 -k 2,3n ${main_dir}/out/${sample}/sv_caller_results${sample}.uwstl_s.nontrs | bgzip -f > ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.nontrs.gz
tabix -p vcf ${main_dir}/out/${sample}/sv_caller_results${sample}.uwstl_s.nontrs.gz
awk '($1=="BND"){print $2"\t"$3"\t"$5"\t"$4"\t"$2":"$3":"$5":"$1}' ${main_dir}/out/${sample}/sv_caller_results/${sample}.chromoseq.txt > ${main_dir}/out/${sample}/sv_caller_results/${sample}.uwstl_s.trs   
         
echo "# prepare the SV vcf files (manta + ichnorcnv)" && date 
gunzip -f ${main_dir}/out/${sample}/sv_caller_results/${sample}.svs_annotated.vcf.gz 
cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.svs_annotated.vcf ${main_dir}/out/${sample}/sv_caller_results/${sample}.manta.vcf         
python ${main_dir}/software/CYTO-SV-ML/Pipeline_script/ichnorcnv_tf.py ${main_dir}/out/${sample}/sv_caller_results/${sample}.segs.txt ${main_dir}/out/${sample}/sv_caller_results/${sample}.ichnorcnv.vcf
