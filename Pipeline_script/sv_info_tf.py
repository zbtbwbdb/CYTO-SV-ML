import sys,os,re
in_vcf=open(sys.argv[1],'r')
out_vcf=open(str(sys.argv[1])+'.tf','w')

def sv_info_tf(line):
    sv_dict={}
    item=line.strip().split('\t')
    sv_chr=item[0]
    sv_chr2=sv_chr
    sv_bp_st=item[1]
    ci_st='-160,160'
    ci_end='-160,160'
    bnd_depth=0
    mate_bnd_depth=0
    
    #info column extraction
    info_list=item[7]
    info=info_list.split(';')      
    for inf in info:
        if inf.startswith('END'):
          sv_bp_end=inf.split('=')[1]
        if inf.startswith('CHR2'):
          sv_chr2=inf.split('=')[1]  
        if inf.startswith('CIPOS'):
          ci_st=inf.split('=')[1] 
        if inf.startswith('CIEND'):
          ci_end=inf.split('=')[1]  
        if inf.startswith('BND_DEPTH'):
          bnd_depth=inf.split('=')[1]  
        if inf.startswith('MATE_BND_DEPTH'):
          mate_bnd_depth=inf.split('=')[1]  
        
    # geno column extraction
    geno_list=item[9]
    geno=info_list.split(':')      
    geno_index=info_list.split(':')      
    gt==geno[geno_index.index('GT')]
    rd==geno[geno_index.index('RD')] 
    rd_pr==geno[geno_index.index('RD')]  
    rd_sp==geno[geno_index.index('RD')] 
    
    sv_id=str(sv_chr)+':'+str(sv_bp_st)+':'+str(sv_bp_end)+':'+str(sv_chr2)          
    line=sv_id+'\t'+str(sv_chr)+'\t'+str(sv_bp_st)+'\t'+str(sv_bp_end)+'\t'+str(sv_chr2)+'\t'+str(gn)+'\t'+str(bnd_depth)+'\t'+str(mate_bnd_depth)+'\t'+str(rd_pr.split(',')[0])+'\t'+str(rd_pr.split(',')[1])+'\t'+str(rd_sp.split(',')[0])+'\t'+str(rd_sp.split(',')[1])+'\t'+str(ci_st.split(',')[0])+'\t'+str(ci_st.split(',')[1])+'\t'+str(ci_end.split(',')[0])+'\t'+str(ci_end.split(',')[1])+'\n'
    return line  
  
  
for line in in_vcf:
    item=line.strip().split('\t')    
    if line.startswith('##'):
        continue
    elif item[0].startswith('#CHROM'):      
        out_vcf.write('sv_id\tsv_chr\tsv_start_bp\tsv_end_bp\tsv_chr2\tsv_type\gn\tbnd_depth\tmate_bnd_depth\trd_pr_ref\trd_pr_alt\trd_sp_ref\trd_sp_alt\tci_st_l\tci_st_r\tci_end_l\tci_end_r\n')
    else:
        out_vcf.write(sv_info_tf(line))        
    
