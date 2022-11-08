################################################################################################################
# TRS SV simplified database format
# ---example----------------------------------------------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_chr2  sv_type  sv_id                                       sv_info
# chr1    100000       1000000    chr21    BND      chr1:100000:1000000:chr21:BND:MANTA*****      ...
# --------------------------------------------------------------------------------------------------------------
# no-TRS SV simplified data format
# ---example----------------------------------------------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_type  sv_id                                                 sv_info
# chr1    10000        1000000    DEL      chr1:10000:1000000:DEL:DELLY*****                       ...
# --------------------------------------------------------------------------------------------------------------
################################################################################################################

import sys,os,re
in_vcf=open(sys.argv[1],'r')
trs_out_vcf=open(str(sys.argv[1])+'.trs','w')
non_trs_out_vcf=open(str(sys.argv[1])+'.nontrs','w')

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
    line=str(item[0])+'\t'+str(item[1])+'\t'+str(END)+'\t'+str(sv_chr2)+'\t'+str(item[0])+':'+str(item[1])+':'+str(END)+':'+str(sv_chr2)+':TRA:'+str(item[2])+'\tTRA\t'+str(item[-1])+'\n'
    return line  
  
def nontrs_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    info_list=item[7]
    info=info_list.split(';')  
    for inf in info:
        if inf.startswith('END'):
            END=inf.split('=')[1]
        if inf.startswith('SVTYPE'): 
            sv_type=inf.split('=')[1]
    # correct sv info        
    if 'SVTYPE' not in info_list:
        if 'DEL|DUP|INV|INS|TRA|BND' in item[2]:
            sv_type=re.findall('DEL|DUP|INV|INS|TRA|BND',item[2])[0]
        elif 'DEL|DUP|INV|INS|TRA|BND' in item[4]:
            sv_type=re.findall('DEL|DUP|INV|INS|TRA|BND',item[4])[0]        
    line=str(item[0])+'\t'+str(item[1])+'\t'+str(END)+'\t'+str(sv_type)+'\t'+str(item[0])+':'+str(item[1])+':'+str(END)+':'+str(item[3])+':'+str(item[2])+'\t'+str(item[0])+'\t'+str(item[-1])+'\n'
    return line  
  
for line in in_vcf:
    item=line.strip().split('\t') 
    if line.startswith('#'):
        continue
#     elif item[0].startswith('#CHROM'):
#         trs_out_vcf.write('sv_chr\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type\tsv_id\tsv_info\n')
#         non_trs_out_vcf.write('sv_chr\tsv_start_bp\tsv_end_bp\tsv_type\tsv_id\tsv_info\n')
    elif 'BND' in item[4] or 'SVTYPE=BND' in item[7] or 'TRA' in item[4] or 'SVTYPE=TRA' in item[7]:
        trs_out_vcf.write(trs_sv_sim(line))        
    else:
        non_trs_out_vcf.write(nontrs_sv_sim(line))         
        
        