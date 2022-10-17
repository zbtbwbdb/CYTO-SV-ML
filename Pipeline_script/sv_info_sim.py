################################################################################################################
# TRS SV simplified data output format
# ---example-----------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_chr2  sv_type  sv_id
# chr1    100000       1000000    chr21    BND      chr1:100000:1000000:chr21:BND:MANTA***** 
# ---------------------------------------------------------------------------
# no-TRS SV simplified data output format
# ---example-----------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_type  sv_id
# chr1    10000        1000000    DEL      chr1:10000:1000000:DEL:DELLY***** 
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
          END=info.split('=')[1]
        if inf.startswith('CHR2'):
          sv_chr2=info.split('=')[1]          
    line=str(item[0])+'\t'+str(item[1])+'\t'+str(END)+'\t'+str(sv_chr2)+'\tBND\t'+str(item[0])+':'+str(item[1])+':'+str(END)+':'+str(sv_chr2)+':BND:'+str(item[2])
    return line  
  
def notrs_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    info_list=item[7]
    info=info_list.split(';')  
    for inf in info:
        if inf.startswith('END'):
          END=info.split('=')[1]
    line=str(item[0])+'\t'+str(item[1])+'\t'+str(END)+'\t'+str(item[3])+'\t'+str(item[0])+':'+str(item[1])+':'+str(END)+':'+str(item[3])+':'+str(item[2])
    return line  
  
for line in in_vcf:
    if not line.startswith('##'):
        continue
    elif item[0].startswith('#'):      
        trs_out_vcf.write(str(sv_chr)+'\t'+str(sv_start_bp)+'\t'+str(sv_end_bp)+'\t'+str(sv_chr2)+'\t'+str(sv_type)+'\t'+str(sv_id))
        non_trs_out_vcf.write(str(sv_chr)+'\t'+str(sv_start_bp)+'\t'+str(sv_end_bp)+'\t'+str(sv_type)+'\t'+str(sv_id))
    elif item[4]=='BND':
        trs_out_vcf.write(trs_sv_sim(line))        
     else:
        non_trs_out_vcf.write(notrs_sv_sim(line))         
        
        
  
