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
        cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.*.vcf  ${main_dir}/out/${sample}/vcf_out/
        gunzip -f  ${main_dir}/out/${sample}/${sample}.svs_annotated.vcf.gz | cat > ${main_dir}/out/${sample}/vcf_out/${sample}.manta.vcf
        python cnv_tf.py ${main_dir}/out/${sample}/${sample}.segs.txt > ${main_dir}/out/${sample}/vcf_out/${sample}.ichnorcnv.vcf
        
        for file in ${main_dir}/out/${sample}/vcf_out/${sample}.*.vcf
            do
                python sv_vcf_tf.py ${file} > ${file}.cr
                svtyper --max_reads 100000 -i ${file} -B ${main_dir}/in/${sample}/${sample}.bam > ${main_dir}/out/${sample}/vcf_out/${file}.svtyper
            done

        # SV sequence complexity run
        for file in ${main_dir}/out/${sample}/vcf_out/${sample}.*.vcf.cr 
            do  
                # make bed file for SV breakpoints
                awk '($1!~"#"){print $0}' ${file} | sed 's% %\t%g' > ${file}.bed
                awk 'FNR==NR{a[$1];b[$1]=$2;next}{c=b[$1]-150 ; if (($2>=150)&&($2<=c)) {$2=$2-150; $3=$2+150; print $0} else if ($2>c) {$2=c-150;$3=c+150; print $0} else if ($2<150){$2=1;$3=300; print $0}}' ${main_dir}/reference/hg38_chromosome_size.txt ${file}.bed | sed 's% %\t%g' > ${file}.bed.bpst
                awk 'FNR==NR{a[$1];b[$1]=$2;next}{$1=$4; $4=$1; c=b[$1]-150 ; if ($3>=c) {$3=c+150;$2=c-150; print $0} else {$2=$3-150; $3=$3+150;  print $0}}' ${main_dir}/reference/hg38_chromosome_size.txt ${file}.bed | sed 's% %\t%g' > ${file}.bed.bpend
                bedtools getfasta -fi ${main_dir}/reference/hg38/hs38.fasta -bed ${file}.bed.bpst -fo ${file}.bed.bpst.fa.out
                bedtools getfasta -fi ${main_dir}/reference/hg38/hs38.fasta -bed ${file}.bed.bpend -fo ${file}.bed.bpend.fa.out

                export PATH=${main_dir}/software/SeqComplex:$PATH
                cd ${main_dir}/software/SeqComplex
                # SV start breakpoint complexity
                perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${file}.bed.bpst.fa.out
                kz --fasta < ${file}.bed.bpst.fa.out > ${file}.bed.bpst.fa.out.kz
                # SV end breakpoint complexity           
                perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${file}.bed.bpend.fa.out
                kz --fasta < ${file}.bed.bpend.fa.out > ${file}.bed.bpend.fa.out.kz            
            done

        # SV calling consoldation
        ls ${main_dir}/out/${sample}/vcf_out/${sample}.*.vcf.svtyper > ${main_dir}/out/${sample}/vcf_out/${sample}.list
        SURVIVOR merge ${main_dir}/out/${sample}/vcf_out/${sample}.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all 

        # SV vcf simplified transformation
        python sv_info_sim.py ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all 
        
        # SV database label
        for SV_database_name in gnomad_qc gnomad_ps 1000g cytoatlas cosmic donor_g
            do
                python sv_database_mapping.py -i ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.nobnd -t ${main_dir}/SV_database/${SV_database_name}.gz -d 1000 -p 0.5 -o ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf_${SV_database_name} 
                python sv_bnd_database_mapping.py ${main_dir}/SV_database/${SV_database_name}.gz ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.bnd  ${SV_database_name} 
                
                # SV info transformation
                python sv_mapping_tf.py ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf_${SV_database_name} 
            done
            
        # SV info consolidate
        python sv_info_tf.py ${main_dir}/out/${sample}/vcf_out/${sample}
        
    done

# combine the sample SV into cohort dataset
sample_all="cohort_name" # please create your own cohort name here
cat ${main_dir}/out/*.sv.all.tf2 > ${main_dir}/out/${sample_all}.sv.all.tf_all

# SV AutoML run
python AutoML.py ${main_dir}/out/${sample_all}.sv.all.tf_all
# Demo: python AutoML.py example/input.csv
