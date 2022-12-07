import sys,os,re
in_vcf=open(sys.argv[1],'r')
trs_out_vcf=open(str(sys.argv[1])+'.trs_tf','w')
nontrs_out_vcf=open(str(sys.argv[1])+'.nontrs_tf','w')

def vcf_bnd_tf(line):
    info_dict={}
    info_list2=[] 
    item=line.strip().split('\t')

    # add CI info
    info=item[7]    
    if not bool(re.search('CIPOS', info)):
        info=str(info)+";CIPOS=-160,160"
    if not bool(re.search('CIEND', info)):  
        info=str(info)+";CIEND=-160,160"  
        
    #info column extraction        
    info_list=info.split(';')             
    i=0
    for inf in info_list:
        if '=' in inf:
            info_dict[inf.split('=')[0]]=inf.split('=')[1]
        else:
            info_list2.append(str(inf))
            i=1
    info_dict.update({'SVTYPE':'BND'})
    #del info_dict["CSQ"]              
    info_dict.pop("CSQ",None)
    
    # geno column extraction
    geno_list=item[9]
    geno=geno_list.split(':')      
    geno_index=item[8].split(':') 
#    print(geno_index)
    gt=geno[geno_index.index('GT')]
    pr=geno[geno_index.index('PR')]
    rd_pr=geno[geno_index.index('PR')].split(',')[1]
    sr="0,0"
    rd_sp=0
    if 'SR' in item[8]:
        sr=geno[geno_index.index('SR')]
        rd_sp=geno[geno_index.index('SR')].split(',')[1] 
    su=int(rd_pr)+int(rd_sp)
    item[8]="GT:PR:SU:PE:SR"
    item[9]=geno[0]+":"+str(pr)+":"+str(su)+":"+str(rd_pr)+":"+str(sr)
    
    # modify info dict and list
    info_dict.update({'SU':su,'PE':rd_pr,'SP':rd_sp})
    info_list2=[]
    for key,value in info_dict.items():
        info_list2.append(str(key)+'='+str(value))
    item[7]=';'.join(str(inf) for inf in info_list2)  

    # print out line   
    line='\t'.join(str(itm) for itm in item)+'\n'
    return line

def vcf_nonbnd_tf(line):
    info_dict={}
    info_list2=[]   
    item=line.strip().split('\t')
    if item[0].startswith('#'):
        return line
    else:
        info=item[7]
        
        # add CI info
        if not bool(re.search('CIPOS', info)):
            info=str(info)+";CIPOS=-160,160"
        if not bool(re.search('CIEND', info)):  
            info=str(info)+";CIEND=-160,160"  
        
        # remove CSQ annotation
        info_list=info.split(';') 
        i=0
        for inf in info_list:
            if '=' in inf:
                info_dict[inf.split('=')[0]]=inf.split('=')[1]
            else:
                info_list2.append(str(inf))
                i=1
        info_dict.pop("CSQ",None)                        
            
        # modify info dict and list
        for key,value in info_dict.items():
            info_list2.append(str(key)+'='+str(value))
        item[7]=';'.join(str(inf) for inf in info_list2)  
        
        line='\t'.join(str(itm) for itm in item)+'\n'
        return line    
  
for line in in_vcf:
    item=line.strip().split('\t')    
    if line.startswith('##FORMAT=<ID=PR'):
        trs_out_vcf.write('##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">\n##FORMAT=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">\n##FORMAT=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">\n') 
        nontrs_out_vcf.write('##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">\n##FORMAT=<ID=SU,Number=.,Type=Integer,Description="Number of pieces of evidence supporting the variant across all samples">\n##FORMAT=<ID=PE,Number=.,Type=Integer,Description="Number of paired-end reads supporting the variant across all samples">\n')         
    elif line.startswith('#'):
        trs_out_vcf.write(line)    
        nontrs_out_vcf.write(line)          
    elif 'SVTYPE=BND' in item[7] or 'MantaBND' in item[7]:
        trs_out_vcf.write(vcf_bnd_tf(line))        
    else:
        nontrs_out_vcf.write(vcf_nonbnd_tf(line))     
