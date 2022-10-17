import sys,os,re
in_vcf=open(sys.argv[1],'r')
out_vcf=open(sys.argv[2],'w')
sm_id=(str(sys.argv[1]).split('\')[-1]).split('.')[0]
                              
def ichnorcnv_tf(line):
    sv_dict={}
    item=line.strip().split('\t')
    if item[0].startswith('ID'):
        line='#CHROM\tPOS\tID\tREF\ALT\QUAL\tFILTER\tINFO\tFORMAT\t'+sm_id+'\n'
        return line
    else:
        sv_type='NA'                    
        genotype='0/1'  
        sv_end=item[3]                               
        if item[11]=='GAIN':
            sv_type='DUP'
            sv_len=int(sv_end)-int(item[2])                              
        elif item[11]=='HETD': 
            sv_type='DEL'  
            sv_len=int(item[2])-int(sv_end)                              
        elif item[11]=='HOMO': 
            sv_type='DEL' 
            sv_len=int(item[2])-int(sv_end)                             
            genotype='1/1' 
        sv_id=':'.join([item[0],item[1],item[2]]) 
        sv_format= 'GT:CN'         
        sv_cnv=item[10]
        log_RD=item[9]
        sv_info= 'END='+str(end)+';SVTYPE='+sv_type+';SVLEN='+sv_len+';IMPRECISE;natorRD='+log_RD                           
        line='\t'.join([item[1],item[2],sv_id,'.',sv_type,'.\t.',sv_info,sv_format,genotype+':'+str(sv_cnv)])+'\n'
        return line  

for line in in_vcf:
    out_vcf.write(ichnorcnv_tf(line))
