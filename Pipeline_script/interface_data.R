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

option_list<-list(make_option(c("-i", "--input_file"),type="character",help="input file", default=NA,metavar="filename"), make_option(c("-o", "--output_file"),type="character",help="output file", default=NA,metavar="filename"))
opt_parser<-OptionParser(option_list=option_list)
opt=parse_args(opt_parser)
input<-as.character(opt$input_file)
opc<-as.character(opt$output_file)

# read in trs and nontrs data
X_trs=read.table(paste(inFile,'trs',sep=''),sep='\t',index=False,header=True)
X_nontrs=read.table(paste(inFile,'nontrs',sep=''),sep='\t',index=False,header=True)
cyto_band_dict=read.table(cyto_dict_file,sep='\t',index=False,header=True)
# combine trs and nontrs data
X_all=rbind(X_trs,X_nontrs)

# save the SV data
save(X_all,cyto_band_dict, file = paste(inFile+'_all',sep=''))
