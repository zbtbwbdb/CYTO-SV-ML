library(optparse)
library(cluster)
library(reghelper)
library(bios2mds)
library(Hmisc)
library(stringr)
library(GOSemSim)
library(org.Hs.eg.db)
library('biomaRt')
################################################################################################################################################################################

option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"), make_option(c("-c", "--cyto_band_file"),type="character",help="cyto band file", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
input<-as.character(opt$input_file)
output<-as.character(opt$output_file)
cyto_band_file<-as.character(opt$cyto_band_file)

# read in trs and nontrs data
X_all=read.table(paste(inFile,'.all_pred',sep=''),sep='\t',index=False,header=True)
cyto_band_dict=read.table(cyto_band_file,sep='\t',index=False,header=True)

# save the SV data
save(X_all,cyto_band_dict, file = output)
