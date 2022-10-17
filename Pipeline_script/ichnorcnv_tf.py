import sys,os,re
in_vcf=open(sys.argv[1],'r')
out_vcf=open(sys.argv[2],'w')
sm_id=(str(sys.argv[1]).split('/')[-1]).split('.')[0]
                              
def ichnorcnv_tf(line):
    sv_dict={}
    item=line.strip().split('\t')
    if item[0].startswith('ID'):
        line='#CHROM\tPOS\tID\tREF\ALT\QUAL\tFILTER\tINFO\tFORMAT\t'+sm_id+'\n'
        return line
    else:                  
        genotype='0/1'  
        sv_end=item[3]                             
        if item[11]=='GAIN' or item[11]=='AMP':
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
        sv_info= 'END='+str(sv_end)+';SVTYPE='+sv_type+';SVLEN='+str(sv_len)+';IMPRECISE;natorRD='+str(log_RD)
        line='\t'.join([str(w) for w in item[1],item[2],sv_id,'.',sv_type,'.\t.',sv_info,sv_format,genotype+':'+str(sv_cnv)])+'\n'
        return line

out_vcf.write('##fileformat=VCFv4.1\n##fileDate=20221013\n##reference=1000GenomesPhase3_decoy-GRCh37\n##source=CNVnator\n##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n##INFO=<ID=natorRD,Number=1,Type=Float,Description="Normalized RD">\n##INFO=<ID=natorP1,Number=1,Type=Float,Description="e-val by t-test">\n##INFO=<ID=natorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">\n##INFO=<ID=natorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">\n##INFO=<ID=natorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">\n##INFO=<ID=natorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">\n##INFO=<ID=natorPE,Number=1,Type=Integer,Description="Number of paired-ends support the event">\n##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">\n##ALT=<ID=DEL,Description="Deletion">\n##ALT=<ID=DUP,Description="Duplication">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">\n##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-ends that support the event">\n') 
for line in in_vcf:    
    item=line.strip().split('\t')
    if item[11]!='NEUT':                         
    	out_vcf.write(ichnorcnv_tf(line))
