import sys,os,re, collections
in_vcf=open(sys.argv[1],'r')
out_vcf=open(str(sys.argv[1])+'.sv_info','w')

def sv_info_tf(line):
    info_dict={}
    item=line.strip().split('\t')
    info_key_list=[]
    
    #info column extraction
    info=item[7]
    info_list=info.split(';') 
    for inf in info_list:
        if '=' in inf:
            info_dict[inf.split('=')[0]]=inf.split('=')[1]
            info_key_list.append(inf.split('=')[0])            
        else:
            info_dict[inf]=inf
    
    # geno column extraction
    geno_index=item[8].split(':')     
    geno_list=item[9].split(':')      
    for i in range(0,len(geno_index)):
        info_dict[geno_index[i]]=geno_list[i]
        
    # correct sv info        
    if 'SVTYPE' not in info_list:
        if 'DEL|DUP|INV|INS|TRA|BND' in item[2]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[2])[0]
        elif 'DEL|DUP|INV|INS|TRA|BND' in item[4]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[4])[0]
            
    alt_info=re.sub('\[|\]',':',item[4])        
    if 'CHR2' not in info_key_list:
        if info_dict['SVTYPE']=='BND':
            for ai in alt_info.split(':'):
                if re.findall('chr',ai):
                    info_dict['CHR2']=ai
        else:
            info_dict['CHR2']=item[0]   
            
    if 'END' not in info_key_list:
        if info_dict['SVTYPE']=='BND':
            for ai in alt_info.split(':'):
                if ai.isdigit():
                    info_dict['END']=ai
        else:
            info_dict['END']=item[1]
                
    # print out dict into line 
    key_list=[]
    value_list=[]
    info_dict = collections.OrderedDict(sorted(info_dict.items()))
    for key,value in info_dict.items():
        key_list.append(str(key))
        value_list.append(str(value))  
    
    if n==0:
        line='\t'.join(str(ch) for ch in col_head)+'\t'+'\t'.join(str(key) for key in key_list)+'\n'+'\t'.join(str(itm) for itm in item[0:6])+'\t'+'\t'.join(str(value) for value in value_list)+'\n'
    else:
        line='\t'.join(str(itm) for itm in item[0:6])+'\t'+'\t'.join(str(value) for value in value_list)+'\n'        
    return line
  
for line in in_vcf:
    item=line.strip().split('\t')
    n=0
    if line.startswith('#CHROM'):
        if '\t' not in line:
            item=line.strip().split(' ')        
        col_head=item[0:6]
        col_head[0]=re.sub("#", "", col_head[0])
        continue     
    elif line.startswith('#'):
        continue           
    else:
        out_vcf.write(sv_info_tf(line))   
        n=1
