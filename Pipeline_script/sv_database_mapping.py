################################################################################################################
# no-TRS SV simplified data format
# ${main_dir}/out/${sample}/vcf_out/${sample}.sv.all.tf.nobnd
# ---example-----------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_type  sv_id
# chr1    10000        1000000    DEL      chr1:10000:1000000:DEL:DellyDEL***** 
# ----------------------------------------------------------------------------
# no-TRS SV simplified database format
# ${main_dir}/SV_database/${SV_database_name}.gz 
# ---example------------------------------------------------------------------
# sv_chr  sv_start_bp  sv_end_bp  sv_type  database_AF
# chr1    10000        1000000    DEL      0.001
# ----------------------------------------------------------------------------
#################################################################################################################

import sys,getopt,os,commands,subprocess,SVCNV,SVCNV_sim,SVCNV_set
#parameter setting
wd = sys.path[0]
opts,args = getopt.getopt(sys.argv[1:],"i:d:p:t:o:")
inFile = ""
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
        if not line.startswith('#') and line[:3]=="chr":
# read basic sv info
            item = line.strip().split("\t")
            chr = item[0]    #chrom
            start = min(int(item[1]),int(item[2]))     #start_pos        
            end = max(int(item[1]),int(item[2]))     #end_pos
            svtype=item[3]  
            sv_id=item[4]            
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
#            print(template_database)
# read database sv_list in the range into dictionary
            svcnv_list = commands.getoutput("tabix -f " + template_database +" "+ chr + ":" + str(start1) + "-" + str(end1) +"|grep "+svtype).split("\n")
#             #svcnv_list = subprocess.check_output("tabix -f " + template_database +" "+ chr + ":" + str(start1) + "-" + str(end1) +"|grep "+svtype).split("\n")
#             process = subprocess.Popen(["tabix -f ",template_database," ",chr, ":",str(start1),"-",str(end1),"|grep ",svtype], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#             out, err = process.communicate()
#             svcnv_list = out.split("\n")
            if len(svcnv_list)==1 and svcnv_list[0]=="":
                fp.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+chr+"\t"+svtype+"\t"+sv_id+"\tNot_in_database\n")
                continue
            else:
                sv_dict=[]
                sv_dict_ori=[]    
                #print(svcnv_list)
                for svcnvs in svcnv_list:
#                    result = svcnvs.strip().split("\t") 
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
                
# check the difference of original sv and total databse svs by mapping out the gap sv          
                if svtype=="DEL" or svtype=="DUP": # or svtype=="INV"
                    sm2= SVCNV_set.subtract_by_overlap(sm1,svcnv1,percent)
                    for sm in sm2:
                        subtract_list_len+=sm.length 
                else:
                    sm2 = SVCNV_set.subtract_by_breakpoint(sm1,svcnv1,distance)
                    for sm in sm2:                      
                        subtract_list_len+=svcnv1.start_pos-sm.start_pos

# calculate the sv consensus by 1 - the difference between original sv and total databse svs by summing up the gap sv length 
            if len(sm2)==0 or svcnv1.length==0:
                svcnv_ratio=1-float(len(filter(None,sm2)))
            else:                        
                svcnv_ratio=float(1-float(subtract_list_len)/float(svcnv1.length)) 

#  determine whether the sv diference is passsing the cutoff and save the data file
            if svcnv_ratio<percent:
                fp.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + svcnv1.chr + "\t"+ svcnv1.svcnv_type + "\t" + sv_id + "\tPASS\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio) + "\t" + "|".join(sv_dict_ori) + "\n" )
            if svcnv_ratio>=percent:
                fp.write(svcnv1.chr + "\t" + str(svcnv1.start_pos) + "\t" + str(svcnv1.end_pos) + "\t" + svcnv1.chr + "\t" + svcnv1.svcnv_type  + "\t" + sv_id + "\tFAIL\t" + svcnv1.info + "\t" + "{:.3f}".format(svcnv_ratio) + "\t" + "|".join(sv_dict_ori) + "\n" )
