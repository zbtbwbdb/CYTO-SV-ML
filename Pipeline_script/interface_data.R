install.packages("optparse")
library(optparse)

################################################################################################################################################################################

option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"), make_option(c("-c", "--cyto_band_file"),type="character",help="cyto band file", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
input<-as.character(opt$input_file)
output<-as.character(opt$output_file)
cyto_band_file<-as.character(opt$cyto_band_file)

# read in trs and nontrs data
all_sv=read.table(paste(input,'.all_pred',sep=''),sep='\t',header=T)
genome_coor=read.table(cyto_band_file,sep='\t',header=T)

# recode the label
all_sv['label']=as.character(all_sv['label'])
all_sv[all_sv['label']=='0', 'label']='UC'
all_sv[all_sv['label']=='-1', 'label']='TA'
all_sv[all_sv['label']=='1', 'label']='TG'
all_sv[all_sv['label']=='2', 'label']='TS'
all_sv['predict_label']=as.character(all_sv['predict_label'])
all_sv[all_sv['predict_label']=='-1','predict_label']='TA'
all_sv[all_sv['predict_label']=='1','predict_label']='TG'
all_sv[all_sv['predict_label']=='2','predict_label']='TS'

# save the SV data
save(all_sv,genome_coor, file = output)
