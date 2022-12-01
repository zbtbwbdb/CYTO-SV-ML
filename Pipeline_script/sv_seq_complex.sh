#!/bin/bash
main_dir=$1
cyto_dv_ml_dir=$2
sample=$3
size_k=$4

# SV sequence complexity run 

echo "# extract SV coordinate information" && date
python ${cyto_dv_ml_dir}/Pipeline_script/sv_vcf_bed_tf.py ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all
awk 'FNR==NR{a[$1];b[$1]=$2;next}{c=b[$1]-150 ; if (($2>=150)&&($2<=c)) {$2=$2-150; $3=$2+150; print $1"\t"$2"\t"$3"\t"$5"\t"$6} else if ($2>c) {$2=c-150;$3=c+150; print $1"\t"$2"\t"$3"\t"$5"\t"$6} else if ($2<150){$2=1;$3=300; print $1"\t"$2"\t"$3"\t"$5"\t"$6}}' ${cyto_dv_ml_dir}/reference/hg38_chromosome_size.txt ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed | sed 's% %\t%g' > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst
awk 'FNR==NR{a[$1];b[$1]=$2;next}{$1=$4; $4=$1; c=b[$1]-150 ; if ($3>=c) {$3=c+150;$2=c-150; print $1"\t"$2"\t"$3"\t"$5"\t"$6} else {$2=$3-150; $3=$3+150;  print $1"\t"$2"\t"$3"\t"$5"\t"$6}}' ${cyto_dv_ml_dir}/reference/hg38_chromosome_size.txt ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed | sed 's% %\t%g' > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpend


export PATH=${main_dir}/software/SeqComplex:$PATH
cd ${cyto_dv_ml_dir}/software/SeqComplex
for bp in bpst bpend
    do
       echo ${bp}
       echo "# make bed file for SV breakpoints" && date       
       bedtools getfasta -fi ${cyto_dv_ml_dir}/reference/hg38/hs38.fasta -bed ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp} -fo ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out
      
       # remove lowcomplex line with "NNNN"+ 1 "A/T/C/G"
       grep NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out -n > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.cr
       cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp
       for i in $(awk '{print $1}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.cr | cut -d ":" -f 1)
           do
               echo $i
               i1=$((i-1))
               awk -v a="$i1" -v b="$i" '{if ((FNR==a)||(FNR==b)) {print $0"\tlc"} else {print $0}}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc
               cp ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp
           done
       awk '($NF!="lc"){print $0}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc
      # rm ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.tmp   
      
       echo "# run SeqComplex and KZ for SV breakpoints" && date
       perl ${main_dir}/software/SeqComplex/profileComplexSeq.pl ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.lc
       ${main_dir}/mambaforge/envs/py27/bin/kz --fasta < ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz
       awk 'FNR==NR{a[$1]; b[$1]=$0; next}{c=$1":"$2"-"$3; if (c in a) {print $0"\t"b[c]}}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp} > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz.index
       awk 'FNR==NR{a[$1]; b[$1]=$0; next}($6 in a) {print $0"\t"b[$6]}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.complex ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz.index > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.${bp}.fa.out.kz.index_complex
    done
  
echo "# consolidate SV breakpoint complexity information" && date  
awk 'FNR==NR{a[$5];b[$5]=$0;next} ($5 in a){print $0"\t"b[$5]}' ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst.fa.out.kz.index_complex ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpend.fa.out.kz.index_complex  > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.fa.out.kz.index_complex
awk 'FNR==NR{a[$5];b[$5]=$0;next} ($6 in a){print $0"\t"b[$6]}'  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.fa.out.kz.index_complex  ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed | awk '{if (FNR==1) {print "CHROM\tPOS\tEND\tsv_chr2\tsvtype\tID\tsv_bp_end_CHROM\tsv_bp_end_POS\tsv_bp_end_END\tsvtype\tID\tsv_bp_end_id\tsv_bp_end_length\tsv_bp_end_cc0\tsv_bp_end_cc1\tsv_bp_end_id\tsv_bp_end_cc_v1\tsv_bp_end_cc_v2\tsv_bp_end_cc_v3\tsv_bp_end_cc_v4\tsv_bp_end_cc_v5\tsv_bp_end_cc_v6\tsv_bp_end_cc_v7\tsv_bp_end_cc_v8\tsv_bp_end_cc_v9\tsv_bp_end_cc_v10\tsv_bp_end_cc_v11\tsv_bp_end_cc_v12\tsv_bp_end_cc_v13\tsv_bp_end_cc_v14\tsv_bp_end_cc_v15\tsv_bp_end_cc_v16\tsv_bp_end_cc_v17\tsv_bp_end_cc_v18\tsv_bp_end_cc_v19\tsv_bp_end_cc_v20\tsv_bp_end_cc_v21\tsv_bp_end_cc_v22\tsv_bp_end_cc_v23\tsv_bp_end_cc_v24\tsv_bp_st_CHROM\tsv_bp_st_POS\tsv_bp_st_END\tsvtype\tID\tsv_bp_st_id\tsv_bp_st_length\tsv_bp_st_cc0\tsv_bp_st_cc1\tsv_bp_st_id\tsv_bp_st_cc_v1\tsv_bp_st_cc_v2\tsv_bp_st_cc_v3\tsv_bp_st_cc_v4\tsv_bp_st_cc_v5\tsv_bp_st_cc_v6\tsv_bp_st_cc_v7\tsv_bp_st_cc_v8\tsv_bp_st_cc_v9\tsv_bp_st_cc_v10\tsv_bp_st_cc_v11\tsv_bp_st_cc_v12\tsv_bp_st_cc_v13\tsv_bp_st_cc_v14\tsv_bp_st_cc_v15\tsv_bp_st_cc_v16\tsv_bp_st_cc_v17\tsv_bp_st_cc_v18\tsv_bp_st_cc_v19\tsv_bp_st_cc_v20\tsv_bp_st_cc_v21\tsv_bp_st_cc_v22\tsv_bp_st_cc_v23\tsv_bp_st_cc_v24\n"$0} else {print $0}}' > ${main_dir}/out/${sample}/${sample}.${size_k}k.sv.all.bed.bpst_bpend.kz.index_complex    
