import sys,os,re
import multiprocessing as mp
vcftyper=open(sys.argv[1],'r')
con=str(sys.argv[2])
#example python sv_count.py test.vcf vc > test.vcf.gene
try:
    no=int(sys.argv[3])
except:
    no=16
    
def vcf_id_count(line):
    item=line.strip().split('\t')
    var_count=0
    if item[0].startswith("#"):
        line='\t'.join(item)
        return line #"NONE"
    else:    
        ref=item[3].split(',')[0]   
        item[3]=ref[:100]        
        alt=item[4].split(',')[0]
        item[4]=alt[:100]
        item[2]=item[0]+':'+item[1]+':'+item[3]+':'+item[4]          
        l_len=len(item)
        for i in range(9,l_len):
            geno=item[i].split(':')[0]
            if not (geno.startswith('./.') or geno.startswith('0/0') or geno.startswith('.|.') or geno.startswith('0|0') or geno=="0"):
                var_count+=1                   
        line='\t'.join(item)+'\t'+str(var_count)
        return line  

test_list=[]                
#pool=mp.Pool(processes=no)
pool=mp.Pool(processes=6)
# id sim and '0/0' correct
if con=="vc":
    output=[pool.apply(vcf_id_count,args=(line,)) for line in vcftyper]
    for line in output:
        if not line.startswith('NONE'):
            print(line)
