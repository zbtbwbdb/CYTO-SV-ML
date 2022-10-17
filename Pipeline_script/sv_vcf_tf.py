import sys,os,re
in_vcf=open(sys.argv[1],'r')
out_vcf=open(str(sys.argv[1])+'.cr','w')

def vcf_sv_tf(line):
    sv_dict={}
    item=line.strip().split('\t')
    if item[0].startswith('#'):
        return line
    else:
        info_list=item[7]
        #info=info_list.split(';')
        if not bool(re.search('CIPOS', info_list)):
            info_list=str(info_list)+";CIPOS=-160,160"
        if not bool(re.search('CIEND', info_list)):  
            info_list=str(info_list)+";CIEND=-160,160"    
        item[7]=info_list
        line='\t'.join(str(w) for w in item)+'\n'
        return line  
    
for line in in_vcf:
    out_vcf.write(vcf_sv_tf(line))
