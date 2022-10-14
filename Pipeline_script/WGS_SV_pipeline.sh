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

# SV preprocess for each sample
for sample in $(cat ${main_dir}/${sample_id_list})
    do
        # SV Parliment2 run
        docker run --rm --privileged -v ${main_dir}/in/${sample}/:/home/dnanexus/in -v ${main_dir}/out/${sample}/:/home/dnanexus/out docker.io/dongwonlee/parliament2-sing:v0.12 --bam ${sample}.bam -r hg38/hs38.fasta --prefix ${sample} --filter_short_contigs --breakdancer --manta --cnvnator --delly_deletion --delly_insertion --delly_inversion --delly_duplication

        # SV ChromoSeq run
        docker run --rm --privileged  -v ${main_dir}/:/scratch  --entrypoint /bin/sh docker.io/zatawada/docker-basespace_chromoseq_v2:master -c '/usr/bin/java -Dconfig.file=/scratch/software/docker-${scratch}/${chromoseq_dir}/lsf/application.new.conf -jar /opt/cromwell-36.jar run -t wdl -i ${scratch}/${chromoseq_dir}/lsf/inputs.json.tmp ${scratch}/${chromoseq_dir}/workflow_files/Chromoseq.v17.wdl'

        # SV svtyper run
        for file in ${main_dir}/out/${sample}/${sample}.*.vcf 
            do
                svtyper --max_reads 100000 -i ${file} -B ${main_dir}/in/${sample}/${sample}.bam > ${file}.svtyper
            done

        # SV sequence complexity run
        for file in ${main_dir}/out/${sample}/${sample}.*.vcf 
            do  
                # make bed file for SV breakpoints
                awk '($1!~"#"){print $0}' ${main_dir}/out/${sample}/${sample}.*.vcf | sed 's% %\t%g' > ${file}.bed
                awk 'FNR==NR{a[$1];b[$1]=$2;next}{c=b[$1]-150 ; if (($2>=150)&&($2<=c)) {$2=$2-150; $3=$2+150; print $0} else if ($2>c) {$2=c-150;$3=c+150; print $0} else if ($2<150){$2=1;$3=300; print $0}}' ${main_dir}/reference/hg38_chromosome_size.txt ${file}.bed | sed 's% %\t%g' > ${file}.bed.bpst
                awk 'FNR==NR{a[$1];b[$1]=$2;next}{$1=$4; $4=$1; c=b[$1]-150 ; if ($3>=c) {$3=c+150;$2=c-150; print $0} else {$2=$3-150; $3=$3+150;  print $0}}' ${main_dir}/reference/hg38_chromosome_size.txt ${file}.bed | sed 's% %\t%g' > ${file}.bed.bpend
                bedtools getfasta -fi ${main_dir}/reference/hg38/hs38.fasta -bed ${file}.bed.bpst -fo ${file}.bed.bpst.fa.out
                bedtools getfasta -fi ${main_dir}/reference/hg38/hs38.fasta -bed ${file}.bed.bpend -fo ${file}.bed.bpend.fa.out

                export PATH=${main_dir}/software/SeqComplex:$PATH
                cd ${main_dir}/software/SeqComplex
                # SV start breakpoint
                perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${file}.bed.bpst.fa.out
                kz --fasta < ${file}.bed.bpst.fa.out > ${file}.bed.bpst.fa.out.kz
                # SV end breakpoint            
                perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${file}.bed.bpend.fa.out
                kz --fasta < ${file}.bed.bpend.fa.out > ${file}.bed.bpend.fa.out.kz            
            done

        # SV calling consoldation
        ls ${main_dir}/out/${sample}/${sample}.*.vcf.svtyper > ${main_dir}/out/${sample}/${sample}.list
        SURVIVOR merge ${main_dir}/out/${sample}/${sample}.list 1000 1 1 0 0 10  ${main_dir}/${sample}/${sample}.sv.all 


        # SV database label
        python sv_database_mapping.py -i ${main_dir}/out/${sample}/${sample}.sv.all.tf.nobnd -t ${main_dir}/SV_database/${SV_database_name}.gz -d 1000 -p 0.5 -o ${main_dir}/out/${sample}/${sample}.sv.all.tf_${SV_database_name} 
        python sv_bnd_database_mapping.py ${main_dir}/SV_database/${SV_database_name}.gz ${main_dir}/out/${sample}/${sample}.sv.all.tf.bnd  ${SV_database_name} ${nf} 

        # SV info transformation


    done

# combine the sample SV into cohort dataset
sample_all="cohort_name" # define your own cohort name here
cat ${main_dir}/out/*.sv.all.tf2 > ${main_dir}/out/${sample_all}.sv.all.tf_all

# SV AutoML run
python AutoML.py ${main_dir}/out/${sample_all}.sv.all.tf_all
# Demo: python AutoML.py example/input.csv
