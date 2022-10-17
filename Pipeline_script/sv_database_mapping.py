################################################################################################################
# no-TRS SV simplified data format
# ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.nobnd
# ---example-----------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_type  sv_id
# chr1    10000        1000000    DEL      chr1:10000:1000000:DEL:DELLY***** 
# ----------------------------------------------------------------------------
# no-TRS SV simplified database format
# ---example------------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_type  database_AF
# chr1    10000        1000000    DEL      0.001
# ----------------------------------------------------------------------------
#################################################################################################################

import sys,getopt,os,commands,SVCNV,SVCNV_sim,SVCNV_set
#parameter setting
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"i:d:p:t:o:")
inFile = ""
cutoff = 0
percent = 0.7
distance = 1000

for op, value in opts:
	if op == "-i":
	    inFile = value
	if op == "-d":
	    distance = int(value)
	if op == "-p":
	    percent = float(value)
	if op == "-t":
	    template_database = value
	if op == "-o":
	    filter_pass = value

if inFile == "":
	print("-i invalid")
	sys.exit()
      
fp=open(filter_pass,'w')

# read sv_list from analysis file
with open(inFile,'r') as f:
	for line in f:
	    if line.startswith('#') or line[:3]!="chr":
#		print(line.strip())
		continue

# read basic sv info
            item = line.strip().split("\t")
            chr = item[0]    #chrom
            start = min(int(item[1]),int(item[2]))     #start_pos        
            end = max(int(item[1]),int(item[2]))     #end_pos
            svtype=item[3]   
            sub_len=int(start-end)

# prepare mapping sv range
            start1 = min(int(start-(end-start)*percent),int(start-1000))
            if start1<0:
                start1=0
            end1 = max(int(end+(end-start)*percent),int(end+1000))

# initiate sv 
            svcnv1=SVCNV_set.SVCNV(line)
            svcnv_list=[]
            sv_dict_ori=[]

# read database sv_list in the range into dictionary
            #print("tabix -f " + chr + ":" + str(start1) + "-" + str(end1) +"|grep "+svtype)
            svcnv_list = commands.getoutput("tabix -f " + template_database +" "+ chr + ":" + str(start1) + "-" + str(end1) +"|grep "+svtype).split("\n")
            if len(svcnv_list)==1 and svcnv_list[0]=="":
                fp.write(line.strip()+"\tNot_in_database\n")
                continue
            else:
                sv_dict=[]
                sv_dict_ori=[]    
                #print(svcnv_list)
                for svcnvs in svcnv_list:
                    result = svcnvs.strip().split("\t")

# check additional database cutoff # remove this because it's currently not standardized                    
#                    if result[4]<=allele_freq or template_database=="/storage/chen/Pipeline/pipeline/pipeline_restructure/trio_cnv/GRCh37_hg19_variants_2016-05-15.txt":
                    #print(svcnvs)
                    s = SVCNV_set.SVCNV(svcnvs)
                    sv_dict.append(s)

# save sv info from class instance as string into a list
                for smo in sv_dict:
                    sv_dict_ori.append(smo.chr + ":" + str(smo.start_pos) + ":" + str(smo.end_pos) + ":" + str(smo.length) + ":" + smo.svcnv_type + ":" + smo.info) 

# simplify sv_dict from database
                if svtype=="DEL" or svtype=="DUP": 
                    sm1=SVCNV_set.simplify_by_overlap(sv_dict)
                else:
                    sm1=SVCNV_set.simplify_by_breakpoint(sv_dict)

# calculate svcnv exclude ratio
                sm2=[]
                sv_dict_sim=[]
                subtract_list_len=0
                svcnv_ratio=0                
                #print(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + svcnv1.svcnv_type + "\t" + str(svcnv1.length))
                
# check the difference of original sv and total databse svs by mapping out the gap sv          
                if svtype=="DEL" or svtype=="DUP":
                    sm2= SVCNV_set.subtract_by_overlap(sm1,svcnv1,percent)
                    for sm in sm2:
                        subtract_list_len+=sm.length 
                else:
                    sm2 = SVCNV_set.subtract_by_breakpoint(sm1,svcnv1,distance)
                    for sm in sm2:                      
                        subtract_list_len+=svcnv1.start_pos-sm.start_pos
 
# calculate the difference of original sv and total databse svs by summing up the gap sv length  


# determine whether the diference of original sv and total databse svs is passsing the cutoff
                if len(sm2)==0 or svcnv1.length==0:    
                    svcnv_ratio=len(filter(None,sm2))
                else:                        
                    svcnv_ratio=float(1-float(subtract_list_len)/float(svcnv1.length)) 
                #print("subtract_list_len:"+str(subtract_list_len)+"svcnv_ratio:"+str(svcnv_ratio))

# save the information into file
                if svcnv_ratio<percent:
                    fp.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + svcnv1.svcnv_type + "\t" + str(svcnv1.length) + "\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio) + "\tPASS\t"+ "|".join(sv_dict_ori) + "\n" )
                if svcnv_ratio>=percent:
                    fp.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + svcnv1.svcnv_type  + "\t" + str(svcnv1.length)+ "\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio) + "\tFAIL\t"+ "|".join(sv_dict_ori) + "\n" )
