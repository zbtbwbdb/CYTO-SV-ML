import sys,os,re
in_vcf=open(sys.argv[1],'r')
out_vcf=open(str(sys.argv[1])+'.re_id','w')
try:
    id_mode=str(sys.argv[2])
except:
    id_mode='a'

def sv_id_tf(line):
    item=line.strip().split('\t')
    sv_chr2=item[0]
    sv_end=item[1]
    sv_type=re.sub('<|>','',item[4])
    
    #collect info from column extraction
    info=item[7]
    info_list=info.split(';') 
    info_dict={}    
    for inf in info_list:
        if '=' in inf:
            info_dict[inf.split('=')[0]]=inf.split('=')[1]
        else:
            info_dict[inf]=inf+"&no&id"       

    # correct sv info        
    if 'SVTYPE=' not in info_list:
        if 'DEL|DUP|INV|INS|TRA|BND' in item[2]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[2])[0]
        elif 'DEL|DUP|INV|INS|TRA|BND' in item[4]:
            info_dict['SVTYPE']=re.findall('DEL|DUP|INV|INS|TRA|BND',item[4])[0]
            
    alt_info=re.sub('\[|\]',':',item[4])     
    if 'CHR2=' not in info_list:
        if info_dict['SVTYPE']=='BND':
            for ai in alt_info.split(':'):
                if re.findall('chr',ai):
                    info_dict['CHR2']=ai
                    sv_chr2=ai
        else:
            info_dict['CHR2']=item[0]   
            
    if 'END=' not in info_list:
        if info_dict['SVTYPE']=='BND':
            for ai in alt_info.split(':'):
                if ai.isdigit():
                    info_dict['END']=ai
                    sv_end=ai
        else:
            info_dict['END']=item[1]
            
    for inf in info_list:        
        if 'CHR2=' in inf:
            sv_chr2=inf.split('=')[1]
        if 'END=' in inf and 'CI' not in inf:
            sv_end=inf.split('=')[1]   
        if 'SVTYPE=' in inf:
            sv_type=inf.split('=')[1]   
            
    # fix sv_id column extraction
    if id_mode=='f' or id_mode=='sf':
        sv_id=item[0]+':'+item[1] +':'+sv_end+':'+sv_chr2 +':'+sv_type
    elif 'chr' in item[2] or re.findall('DEL|DUP|INV|INS|TRA|BND',item[2]) or re.findall('del|dup|inv|ins|tra|bnd',item[2]) or ':' in item[2]:
        sv_id=item[2]
    else:
        sv_id=item[0]+':'+item[1] +':'+sv_end+':'+sv_chr2 +':'+sv_type
    item[2]=sv_id
    
    # re-compose info columns
    info_new=""
    m=0
    for key,value in info_dict.items():
        if "&no&id" in key and m==0 and key!="CSQ":
            info_new=value
            m=1
        elif "&no&id" in key and m==1 and key!="CSQ":
            info_new=info_new+";"+value
            m=1
        elif "&no&id" not in key and m==0 and key!="CSQ":
            info_new=str(key)+"="+str(value) 
            m=1
        elif "&no&id" not in key and m==1 and key!="CSQ":
            info_new=info_new+";"+str(key)+"="+str(value) 
            m=1   
    print(info_new)
    if ';CSQ=' in item[7] and id_mode!='sf': 
        item[7]=info_new+';CSQ='+info_dict['CSQ']
    
    # print out line
    line='\t'.join(str(w) for w in item)+'\n'        
    return line
  
for line in in_vcf:
    item=line.strip().split('\t')
    n=0
    if line.startswith('#'):
        out_vcf.write(line)        
    else:
        out_vcf.write(sv_id_tf(line))   
