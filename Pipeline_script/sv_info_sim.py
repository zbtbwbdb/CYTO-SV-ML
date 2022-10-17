################################################################################################################
# TRS SV simplified data output format
# ---example-----------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_chr2  sv_type  sv_id
# chr1    100000       1000000    chr21    BND      chr1:100000:1000000:chr21:BND:MantaBND***** 
# ---------------------------------------------------------------------------
# no-TRS SV simplified data output format
# ---example-----------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_type  sv_id
# chr1    10000        1000000    DEL      chr1:10000:1000000:DEL:DellyDEL***** 
# ----------------------------------------------------------------------------
################################################################################################################

import sys,os,re
in_vcf=open(sys.argv[1],'r')
trs_out_vcf=open(str(sys.argv[1])+'.tf.trs','w')
non_trs_out_vcf=open(str(sys.argv[1])+'.tf.notrs','w')

def trs_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    info_list=item[7]
    info=info_list.split(';')  
    for inf in info:
        if inf.startswith('END'):
          END=inf.split('=')[1]
        if inf.startswith('CHR2'):
          sv_chr2=inf.split('=')[1]          
    line=str(item[0])+'\t'+str(item[1])+'\t'+str(END)+'\t'+str(sv_chr2)+'\tBND\t'+str(item[0])+':'+str(item[1])+':'+str(END)+':'+str(sv_chr2)+':BND:'+str(item[2]+'\n')
    return line  
  
def notrs_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    info_list=item[7]
    info=info_list.split(';')  
    for inf in info:
        if inf.startswith('END'):
          END=inf.split('=')[1]
    line=str(item[0])+'\t'+str(item[1])+'\t'+str(END)+'\t'+str(item[3])+'\t'+str(item[0])+':'+str(item[1])+':'+str(END)+':'+str(item[3])+':'+str(item[2]+'\n')
    return line  
  
for line in in_vcf:
    item=line.strip().split('\t')    
    if line.startswith('##'):
        continue
    elif item[0].startswith('#CHROM'):      
        trs_out_vcf.write('sv_chr\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type\tsv_id\n')
        non_trs_out_vcf.write('sv_chr\tsv_start_bp\tsv_end_bp\tsv_type\tsv_id\n')
    elif item[4]=='BND':
        trs_out_vcf.write(trs_sv_sim(line))        
    else:
        non_trs_out_vcf.write(notrs_sv_sim(line))         
        
        
  
