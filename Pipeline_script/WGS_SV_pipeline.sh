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

        # SV vcf preparation
        cp ${main_dir}/out/${sample}/sv_caller_results/${sample}.*.vcf  ${main_dir}/out/${sample}/vcf_out/
        gunzip -f  ${main_dir}/out/${sample}/${sample}.svs_annotated.vcf.gz | cat > ${main_dir}/out/${sample}/vcf_out/${sample}.manta.vcf
        python ichnorcnv_tf.py ${main_dir}/out/${sample}/${sample}.segs.txt ${main_dir}/out/${sample}/vcf_out/${sample}.ichnorcnv.vcf
        
        # SV vcf consoldation
        ls ${main_dir}/out/${sample}/vcf_out/${sample}.*.vcf.svtyper > ${main_dir}/out/${sample}/vcf_out/${sample}.list
        SURVIVOR merge ${main_dir}/out/${sample}/vcf_out/${sample}.list 1000 1 1 0 0 10  ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all 
        
        # SV svtyper run           
        python sv_vcf_tf.py ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all
        svtyper --max_reads 100000 -i ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.cr -B ${main_dir}/in/${sample}/${sample}.bam > ${main_dir}/out/${sample}/vcf_out/${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.svtyper

        # SV sequence complexity run 
        # make bed file for SV breakpoints
        awk '($1!~"#"){print $0}' ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.cr | sed 's% %\t%g' > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed
        awk 'FNR==NR{a[$1];b[$1]=$2;next}{c=b[$1]-150 ; if (($2>=150)&&($2<=c)) {$2=$2-150; $3=$2+150; print $0} else if ($2>c) {$2=c-150;$3=c+150; print $0} else if ($2<150){$2=1;$3=300; print $0}}' ${main_dir}/reference/hg38_chromosome_size.txt ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed | sed 's% %\t%g' > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpst
        awk 'FNR==NR{a[$1];b[$1]=$2;next}{$1=$4; $4=$1; c=b[$1]-150 ; if ($3>=c) {$3=c+150;$2=c-150; print $0} else {$2=$3-150; $3=$3+150;  print $0}}' ${main_dir}/reference/hg38_chromosome_size.txt ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed | sed 's% %\t%g' > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpend
        bedtools getfasta -fi ${main_dir}/reference/hg38/hs38.fasta -bed ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpst -fo ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpst.fa.out
        bedtools getfasta -fi ${main_dir}/reference/hg38/hs38.fasta -bed ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpend -fo ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpend.fa.out

        export PATH=${main_dir}/software/SeqComplex:$PATH
        cd ${main_dir}/software/SeqComplex
        # SV start breakpoint complexity
        perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpst.fa.out
        kz --fasta < ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpst.fa.out > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpst.fa.out.kz
        # SV end breakpoint complexity           
        perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpend.fa.out
        kz --fasta < ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpend.fa.out > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.bed.bpend.fa.out.kz    


        # SV vcf simplified transformation
        python sv_info_sim.py ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all 
        awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $6"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5}}' ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.trs > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.trs_ano
        awk '{if (FNR==1) {print "sv_id\tsv_chr1\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type"} else {print $5"\t"$1"\t"$2"\t"$3"\t"$1"\t"$4}}' ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.notrs > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.notrs_ano
        
        # SV database annotation label
        for SV_database_name in gnomad_qc gnomad_ps 1000g cytoatlas cosmic donor_g
            do
                python sv_database_mapping.py -i ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.notrs -t ${main_dir}/SV_database/${SV_database_name}.nontrs.gz -d 1000 -p 0.5 -o ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.notrs_${SV_database_name} 
                python sv_bnd_database_mapping.py ${main_dir}/SV_database/${SV_database_name}.trs.gz ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.trs ${SV_database_name} 
                
                # SV database annotation consolidation/transformation
                python sv_db_tf.py ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.trs_ano ${SV_database_name} 
                python sv_db_tf.py ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.notrs_ano ${SV_database_name}   
            done
            
        # combine the sv annotation vcf           
        cat ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.trs_ano ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.notrs_ano > ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.all_ano
                
        # combine the sv annotation and complexity and svtyper info
        python sv_info_tf.py ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.all_ano          

        # combine the sample SV into cohort dataset
        sample_all="cohort_name" # please create your own cohort name here
        cat ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.all_ano >> ${main_dir}/out/${sample_all}.sv.all.tf.all_ano
    done

# SV AutoML run
python AutoML.py ${main_dir}/out/${sample_all}.sv.all.tf.all_ano
# Demo: python AutoML.py example/input.csv
