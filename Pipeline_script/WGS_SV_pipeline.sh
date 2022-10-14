# pre-installation
# Parliment2: docker pull dongwonlee/parliament2-sing:v0.12
# ChromoSeq: docker pull zatawada/docker-basespace_chromoseq_v2:master
SURVIVOR: https://github.com/fritzsedlazeck/SURVIVOR
SVTyper: https://github.com/hall-lab/svtyper
SeqComplex: https://github.com/caballero/SeqComplex
Komplexity: https://github.com/eclarke/komplexity
#################################################################################################################
#!/bin/bash

# Parliment2 run
sudo docker run --rm --privileged -v ${main_dir}/in/:/home/dnanexus/in -v ${main_dir}/out/:/home/dnanexus/out  docker.io/dongwonlee/parliament2-sing:v0.12 --bam ${sample}.bam --bai ${sample}.bam.bai -r hg38/hs38.fasta --fai  hg38/hs38.fasta.fai --prefix  ${sample} --filter_short_contigs  --breakdancer --breakseq --manta --cnvnator     --lumpy --delly_deletion --delly_insertion --delly_inversion --delly_duplication --genotype 

# ChromoSeq run
sudo docker run --rm --privileged  -v ${main_dir}/:/scratch  --entrypoint /bin/sh docker-basespace_chromoseq:latest -c '/usr/bin/java -Dconfig.file=/scratch/software/docker-${scratch}/${chromoseq_dir}/lsf/application.new.conf -jar /opt/cromwell-36.jar run -t wdl -i ${scratch}/${chromoseq_dir}/lsf/inputs.json.tmp ${scratch}/${chromoseq_dir}/workflow_files/Chromoseq.v17.wdl'

# svtyper run


# Sequence complexity run


# SV calling consoldation
SURVIVOR merge ${outdir}/${sample}/${sample}.list 1000 1 1 0 0 10  ${outdir}/${sample}/${sample}.SV.all 

# SV database label
python ../SVCNV_database_filter.py -i ${outdir}/${sample}/${sample}.SV.all.tf -t ${outdir}/SV_database/${SV_database} -d 1000 -p 0.5 -o ${outdir}/${sample}/${sample}.SV.all.tf_${SV_database} 
python sv_bnd_info_mapping.py ${outdir}/${sample}/${sample}.SV.all.tf ${outdir}/SV_database/${SV_database} ${SV_database} ${nf} 

