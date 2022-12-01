#!/bin/bash
main_dir=$1
parliment_docker=$2
sample=$3  

echo "# run parliament docker" && date
mkdir ${main_dir}/in/reference/
cp ${main_dir}/in/${sample}.bam ${main_dir}/in/input.bam
cp ${main_dir}/in/${sample}.bam.bai ${main_dir}/in/input.bam.bai
cp ${main_dir}/reference/hg38/hg38.fasta ${main_dir}/in/reference/hg38/hg38.fasta
cp ${main_dir}/reference/hg38/hg38.fasta.fai ${main_dir}/in/reference/hg38/hg38.fasta.fai
sudo docker run --rm --privileged -v ${main_dir}/in/:/home/dnanexus/in -v ${main_dir}/out/${sample}/:/home/dnanexus/out ${parliment_docker} --bam ${sample}.bam  --bai ${sample}.bam.bai -r reference/hg38/hg38.fasta --fai reference/hg38/hg38.fasta.fai --prefix ${sample} --filter_short_contigs --breakdancer --cnvnator --delly_deletion --delly_inversion --delly_duplication
rm -rf ${main_dir}/in/input.bam*
rm -rf ${main_dir}/in/reference/
