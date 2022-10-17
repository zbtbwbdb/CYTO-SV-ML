import sys,os,re
in_vcf=open(sys.argv[1],'r')
out_vcf=open(sys.argv[2],'w')

def vcf_sv_sim(line):
    sv_dict={}
    item=line.strip().split('\t')
    if item[0].startswith('#'):
#        line='\t'.join(item)+'\n'
        return line
    else:
        info_list=item[7]
        info=info_list.split(';')
        for inf in info:
            if "=" in inf:
                sv_dict[inf.split('=')[0]]=inf.split('=')[1]
        id=':'.join([item[0],item[1],sv_dict['END'],sv_dict['SVTYPE'],item[2]]) 
        line='\t'.join([item[0],item[1],sv_dict['END'],sv_dict['SVTYPE'],id])+'\n'
    #        line='\t'.join([item[0],item[1],sv_dict['END'],sv_dict['SVTYPE'],item[7]])+'\n'
        return line  

for line in in_vcf:
#    if not line.startswith('##'):
    out_vcf.write(vcf_sv_sim(line))
