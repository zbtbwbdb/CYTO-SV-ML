#!/bin/bash
main_dir=$1
sample=$2   
gender=$3

echo "# set up the configure file for chromoseq" && date
sed 's%/scratch/out/XXXXXX%/scratch/out/${sample}/sv_caller_results/%g' ${main_dir}/software/docker-basespace_chromoseq/lsf/inputs.json | sed "s%XXXXXX%${sample}%g"  | sed "s%Male%${gender}%g" > ${main_dir}/software/docker-basespace_chromoseq/lsf/inputs.json.tmp

echo "# run chromoseq docker" && date
sudo docker run --rm --privileged  -v ${main_dir}/:/scratch  --entrypoint /bin/sh zatawada/docker-basespace_chromoseq_v2:master -c '/usr/bin/java -Dconfig.file=/scratch/software/docker-basespace_chromoseq/lsf/application.new.conf -jar /opt/cromwell-36.jar run -t wdl -i /scratch/software/docker-basespace_chromoseq/lsf/inputs.json.tmp /scratch/software/docker-basespace_chromoseq/workflow_files/Chromoseq.v17.wdl'     
