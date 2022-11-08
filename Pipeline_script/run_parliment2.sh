#!/bin/bash
main_dir=$1
sample=$2   

echo "# run parliament docker" && date
sudo docker run -u root --rm --privileged -v ${main_dir}/in/:/home/dnanexus/in -v ${main_dir}/out/${sample}/:/home/dnanexus/out docker.io/dnanexus/parliament2:latest --bam ${sample}.bam  --bai ${sample}.bam.bai -r reference/hg38/hs38.fasta --fai reference/hg38/hs38.fasta.fai --prefix ${sample} --filter_short_contigs --breakdancer --cnvnator --delly_deletion --delly_inversion --delly_duplication
#sudo docker run -u "$(id -u):$(id -g)" --rm --privileged -v ${main_dir}/in/:/home/dnanexus/in -v ${main_dir}/out/${sample}/:/home/dnanexus/out docker.io/dnanexus/parliament2:latest --bam ${sample}.bam  --bai ${sample}.bam.bai -r reference/hg38/hs38.fasta --fai reference/hg38/hs38.fasta.fai --prefix ${sample} --filter_short_contigs --breakdancer --cnvnator --delly_deletion --delly_inversion --delly_duplication
#sudo docker run --rm --privileged -v ${main_dir}/in/:/home/dnanexus/in -v ${main_dir}/out/${sample}/:/home/dnanexus/out docker.io/dongwonlee/parliament2-sing:v0.12 --bam ${sample}.bam  --bai ${sample}.bam.bai -r reference/hg38/hs38.fasta --fai reference/hg38/hs38.fasta.fai --prefix ${sample} --filter_short_contigs --breakdancer --cnvnator --delly_deletion --delly_inversion --delly_duplication